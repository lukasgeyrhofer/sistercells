#!/usr/bin/env python

import numpy as np
import pandas as pd

class SisterCellData(object):
    def __init__(self,**kwargs):
        self.__infiles    = kwargs.get('infiles',[])
        self.__debugmode  = kwargs.get('debugmode',False) # return only single trajectory in iteration
        self.__sisterdata = kwargs.get('sisterdata',True)
        
        # lists to store data internally
        self.__data = list()
        self.__dataorigin = list()
        self.__keylist = list()
                
        # load first sheet of each Excel-File, fill internal data structure
        for filename in self.__infiles:
            try:
                tmpdata = pd.read_excel(filename)
            except:
                continue
            self.__data.append(tmpdata)
            self.__dataorigin.append(filename)
            for k in tmpdata.keys():
                if not str(k) in self.__keylist:
                    self.__keylist.append(str(k))

        
        # there's no point in not having data ...
        # ... or something went wrong. rather stop here
        if not len(self) > 0:
            raise IOError('no data loaded')


    def LMSQ(self,x,y,cov = True):
        """
        Least Mean Squares estimator
        for a linear interpolation between x,y: y ~ a + b x
        Returns also covariance matrix for the estimated parameters a,b if 'cov' is True
        
        Algorithm:
        A = ( 1  ... 1  )
            ( x1 ... xn )
            
        y = ( y1 ... yn ).T
        
        p = (a b).T
        
        E[p]     = inv(A.T * A) * A.T * y
        Cov[p,p] = sigma2 * inv(A.T * A)
        sigma2   = E[ ( y - A*E[p] )^2 ]
        """
        
        n   = len(x)
        sx  = np.sum(x)
        sy  = np.sum(y)
        sxx = np.dot(x,x)
        sxy = np.dot(x,y)
        syy = np.dot(y,y)
        
        # estimate parameters
        denom     = (n*sxx-sx*sx)
        b         = (n*sxy - sx*sy)/denom
        a         = (sy-b*sx)/n
        estimate  = np.array([a,b],dtype=np.float)

        if cov:
            # estimate covariance matrix of estimated parameters
            sigma2    = syy + n*a*a + b*b*sxx + 2*a*b*sx - 2*a*sy - 2*b*sxy # variance of deviations from linear line
            covmatrix = sigma2 / denom * np.array([[sxx,-sx],[-sx,n]],dtype=np.float)

            return estimate,covmatrix
        else:
            return estimate



    def otsu(self,x):
        """
        Otsu's method
        First described in IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS (1979)
        Usually used to binarize photos into black/white
        Here modified to use bins with single entries for each measurement
        """
        sx = np.sort(x)
        lx = len(x)

        w   = np.arange(lx,dtype=np.float)/lx
        m   = np.array([np.mean(sx[:k]) * w[k] if k > 0 else 0 for k in range(lx)])
        mT  = np.mean(sx)
        
        sB  = np.array([(mT * w[k] - m[k])**2/(w[k]*(1.-w[k])) if w[k]*(1.-w[k]) > 0 else 0 for k in range(lx)])
        idx = np.argmax(sB)
        
        return sx[idx]




    def CellDivisionTrajectory(self,dataID, discretize_by = 'length', sisterdata = None, additional_columns = []):
        """
        Transform measured time series data into data for each generation:
        (1) Find cell division events as a large enough drop in the measured value given by 'discretize_by'
        (2) Compute various observables from these generations of cells
        (3) Returns two pandas dataframes for each of the discretized trajectories
        """
        if not discretize_by in self.keylist_stripped:  raise KeyError('key not found')
        # not sure if this is needed
        if sisterdata is None:  sisterdata = self.__sisterdata
        if sisterdata:          keysuffix = ['A','B']    
        else:                   keysuffix = ['']
        
        # return two pandas-dataframes for the two trajectories usually contained in one file as two elements of a list
        # as the two sisters do not need to have the same number of cell divisions, a single dataframe might cause problems with variably lengthed trajectories
        ret = list()
        
        alldiffdata = np.concatenate([np.diff(self.__data[dataID][discretize_by + ks]) for ks in keysuffix])
        
        for ks in keysuffix:
            # use 'discretize_by' as the column name that should serve as measurement that indicates the division
            # in our data, both 'length' and 'cellsize' should be OK.
            # get list of changes between measurements. Usually this distribution is bimodal:
            #  * scattering around small values for normal data-points
            #  * scattering around 'large' negative values for cell divisions
            diffdata  = np.diff(self.__data[dataID][discretize_by + ks])
            # estimate a threshold from the data between the two peaks in the bimodal distribution with Otsu's method, then get indices of these transitions
            index_div = np.where(diffdata < self.otsu(alldiffdata))[0].flatten()
            
            # 'double jump': division spans two (or more) time points
            # todo: need to implement some lines of code that detect if a division spans multiple time points, and corrects the algorithm below properly
            dj = sum(np.where(np.diff(index_div) == 1)[0].flatten())
            
            # timepoint of division is assumed to be the average before and after the drop in signal
            time_div  = np.concatenate([[1.5*self.__data[dataID]['time' + ks][0] - 0.5 * self.__data[dataID]['time' + ks][1]],0.5 * np.array(self.__data[dataID]['time' + ks][index_div + 1]) + 0.5 * np.array(self.__data[dataID]['time'+ks][index_div])])
            
            # store results in dictionary, which can be easier transformed into pandas-dataframe
            # values computed for each cell cycle are
            #  * generation time
            #  * size at birth
            #  * size at division
            #  * Least Mean Squares Estimator for the (exponential) growth rate over the full division cycle
            ret_ks = dict()
            #ret_ks['divisiontimes']           = np.array(time_div)
            ret_ks['generationtime']          = np.diff(time_div)
            ret_ks[discretize_by + '_birth']  = np.concatenate([[self.__data[dataID][discretize_by + ks][0]],np.array(self.__data[dataID][discretize_by + ks][index_div+1])[:-1]])
            ret_ks[discretize_by + '_final']  = np.array(self.__data[dataID][discretize_by + ks][index_div])
            ret_ks['growth_' + discretize_by] = np.concatenate([
                            [self.LMSQ(self.__data[dataID]['time' + ks][:index_div[0]],np.log(self.__data[dataID][discretize_by + ks][:index_div[0]]),cov=False)[1]], # first generation
                            np.array([self.LMSQ(
                                self.__data[dataID]['time' + ks][index_div[i]+1:index_div[i+1]+1],                  # x-values
                                np.log(self.__data[dataID][discretize_by + ks][index_div[i]+1:index_div[i+1]+1]),   # y-values
                                cov=False)[1] for i in range(len(index_div)-1)])])
            
            # if additional data is requested, return values from birth and division
            if len(additional_columns) > 0:
                for ac in additional_columns:
                    if ac in self.keylist_stripped:
                        ret_ks[ac + '_birth'] = np.concatenate([[self.__data[dataID][ac + ks][0]],np.array(self.__data[dataID][ac + ks][index_div+1])[:-1]])
                        ret_ks[ac + '_final'] = np.array(self.__data[dataID][ac + ks][index_div])
            
            
            # we have everything, now make a dataframe
            ret.append(pd.DataFrame(ret_ks))
        return ret


    def AutocorrelationFFT(self, dataID, normalize = False, maxlen = None, enforcelen = False) :
        """
        Compute the autocorrelation of the signal, based on the properties of the
        power spectral density of the signal.
        """
        trajA, trajB = self.CellDivisionTrajectory(dataID)

        acfA = dict()
        for key in trajA.keys():
            ml  = np.min([maxlen,len(trajA[key])])
            xp  = trajA[key][:ml]-np.mean(trajA[key][:ml])
            f   = np.fft.fft(xp)
            p   = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
            pi  = np.fft.ifft(p)
            acf = np.real(pi)[:ml/2]
            if normalize: acf *= 1./np.sum(xp**2)
            if enforcelen and not maxlen is None:
                if len(acf) < maxlen:
                    acf = np.concatenate([acf,np.zeros(maxlen - len(acf))])
            acfA[key] = acf

        for key in trajB.keys():
            ml  = np.min([maxlen,len(trajB[key])])
            xp  = trajB[key][:ml]-np.mean(trajB[key][:ml])
            f   = np.fft.fft(xp)
            p   = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
            pi  = np.fft.ifft(p)
            acf = np.real(pi)[:ml/2]
            if normalize: acf *= 1./np.sum(xp**2)
            if enforcelen and not maxlen is None:
                if len(acf) < maxlen:
                    acf = np.concatenate([acf,np.zeros(maxlen - len(acf))])
            acfB[key] = acf
            
        return acfA,acfB


    def AutocorrelationRestricted(self, dataID, maxlen = 20):
        """
        Compute autocorrelation on same set of data as the lineagecorrelation below
        (only use first 'maxlen' steps of discretized trajectory)
        """
        
        trajA, trajB = self.CellDivisionTrajectory(dataID)
        
        acf_sumAA   = dict()
        acf_sumA    = dict()
        acf_countAA = dict()
        acf_countA  = dict()
        
        for corrkey in trajA.keys():
            if not corrkey in acf_sumAA.keys():
                acf_sumAA  [corrkey] = np.zeros(maxlen, dtype = np.float)
                acf_countAA[corrkey] = np.zeros(maxlen, dtype = np.float)
                acf_sumA   [corrkey] = 0
                acf_countA [corrkey] = 0
            
            mlA = min(maxlen, len(trajA))
            mlB = min(maxlen, len(trajB))

            acf_sumA[corrkey]          += np.sum(trajA[corrkey][:mlA])
            acf_countA[corrkey]        += mlA
            
            acf_sumA[corrkey]          += np.sum(trajB[corrkey][:mlB])
            acf_countA[corrkey]        += mlB
            
            acf_sumAA[corrkey][:mlA]   += trajA[corrkey][0] * trajA[corrkey][:mlA]
            acf_countAA[corrkey][:mlA] += np.ones(mlA)
            
            acf_sumAA[corrkey][:mlB]   += trajB[corrkey][0] * trajB[corrkey][:mlB]
            acf_countAA[corrkey][:mlB] += np.ones(mlB)
            
        return acf_sumAA, acf_countAA, acf_sumA, acf_countA



    def LineageCorrelation(self, dataID, maxlen = 20):
        """
        Compute the correlation between the two lineages A and B of sistercells,
        specifically <X(A,t) X(B,t)> - <X(A,t)> <X(B,t)>,
        where X indicates the data in the discretized trajectory
        """
        trajA,trajB = self.CellDivisionTrajectory(dataID)
        corrmatrix_sumAB = dict()
        corrmatrix_sumA  = dict()
        corrmatrix_sumB  = dict()
        corrmatrix_count = dict()
        
        for corrkey in trajA.keys():
            if not corrkey in corrmatrix_sumAB.keys():
                corrmatrix_sumAB[corrkey] = np.zeros((maxlen,maxlen), dtype = np.float)
                corrmatrix_sumA [corrkey] = np.zeros((maxlen,maxlen), dtype = np.float)
                corrmatrix_sumB [corrkey] = np.zeros((maxlen,maxlen), dtype = np.float)
                corrmatrix_count[corrkey] = np.zeros((maxlen,maxlen), dtype = np.float)
        
            mlA = min(maxlen, len(trajA))
            mlB = min(maxlen, len(trajB))
            
            corrmatrix_sumAB[corrkey][:mlA,:mlB] += np.outer(    trajA[corrkey][:mlA], trajB[corrkey][:mlB] )
            corrmatrix_sumB [corrkey][:mlA,:mlB] += np.repeat( [ trajB[corrkey][:mlB] ], mlA, axis = 0)
            corrmatrix_sumA [corrkey][:mlA,:mlB] += np.repeat( [ trajA[corrkey][:mlA] ], mlB, axis = 0).T
            corrmatrix_count[corrkey][:mlA,:mlB] += np.ones( (mlA, mlB) )
                    
        return corrmatrix_sumAB, corrmatrix_sumA, corrmatrix_sumB, corrmatrix_count
    

    def DivisionTimeDifferences(self,dataID):
        """
        compute the differences in division time points,
        and how this difference progresses in the measurement
        """
        trajA, trajB = self.CellDivisionTrajectory(dataID)
        trajlen = np.min([len(trajA['generationtime']), len(trajB['generationtime'])])

        dt = np.zeros(trajlen)
        for i in range(1,trajlen):
            dt[i] = dt[i-1] + trajA['generationtime'][i] - trajB['generationtime'][i]
        
        return dt



    # control: pair up trajectories randomly
    def RandomTrajectories(self,count = 1, force_length = True):
        for i in range(count):
            tID = np.random.randint(len(self.__data), size = 2)
            ab  = [chr(65+x) for x in np.random.randint(2, size = 2)]
            keys0 = [k for k in self[tID[0]].keys() if k[-1] == ab[0]]
            keys1 = [k for k in self[tID[1]].keys() if k[-1] == ab[1]]
            
            if force_length:
                mlen  = min([len(self[tID[0]]),len(self[tID[1]])])
                newdf = pd.concat([self[tID[0]][keys0][:mlen], self[tID[1]][keys1][:mlen]], axis = 1)
            else:
                newdf = pd.concat([self[tID[0]][keys0], self[tID[1]][keys1]], axis = 1)
                
            newdf.columns = np.concatenate([[k[:-1] + 'A' for k in keys0],[k[:-1] + 'B' for k in keys1]])
            
            yield newdf



    # access single dataframe by its ID
    def __getitem__(self,key):
        return self.__data[key]

    
    # should not allow accessing internal variables in other ways than funneling through this here
    def __getattr__(self,key):
        if key == "filenames":
            return self.__dataorigin
        elif key == "keylist":
            return self.__keylist
        elif key == "keylist_stripped":
            return list(set([s.strip('AB ') for s in self.__keylist]))
        elif key == 'timestep':
            ts = self[0]['timeA'][1] - self[0]['timeA'][0]
            for dataID in range(len(self)):
                dt = np.diff(self.__data[dataID]['timeA'])
                if ts > np.min(dt): ts = np.min(dt)
            return ts

    
    # convenience
    def __len__(self):
        return len(self.__data)


    # data will be processes as loop over the class instance
    # 'debugmode' only returns a single item (the first one)
    def __iter__(self):
        dataIDs = np.arange(len(self),dtype=int)
        if self.__debugmode:
            # yield only first element in debugmode to check analysis on this trajectory
            yield 0,self.__dataorigin[0],self.__data[0]
        else:
            for dataID,origin,data in zip(dataIDs,self.__dataorigin,self.__data):
                yield dataID,origin,data




        
