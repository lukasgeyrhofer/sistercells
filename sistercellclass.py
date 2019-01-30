#!/usr/bin/env python

import numpy as np
import pandas as pd

class SisterCellData(object):
    def __init__(self,**kwargs):
        self.__infiles    = kwargs.get('infiles',[])
        self.__debugmode  = kwargs.get('debugmode',False) # return only single trajectory in iteration
        
        # lists to store data internally
        self.__data = list()
        self.__dataorigin = list()
        self.__keylist = list()

        # load first sheet of each Excel-File, fill internal data structure
        for filename in self.__infiles:
            tmpdata = pd.read_excel(filename)
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
        # least mean squares estimator
        # for a linear interpolation, including covariance matrix for the estimated parameters
        # A = ( 1  ... 1  )
        #     ( x1 ... xn )
        #
        # y = ( y1 ... yn ).T
        #
        # p = (a b).T
        #
        # E[p]     = inv(A.T * A) * A.T * y
        # Cov[p,p] = sigma2 * inv(A.T * A)
        # sigma2   = E[ ( y - A*E[p] )^2 ]
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



    def otsu(self,x,historange = None, bins = None):
        # otsu's method
        # described in IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS (1979)
        # usually used to binarize photos into black/white
        if historange is None:
            if bins is None:    count,binedges = np.histogram(x)
            else:               count,binedges = np.histogram(x,bins=bins)
        else:
            if bins is None:    count,binedges = np.histogram(x,range=historange)
            else:               count,binedges = np.histogram(x,range=historange,bins=bins)
        bincenter = binedges[:-1] + .5 * np.diff(binedges)
        bins      = len(bincenter)

        p   = count/float(sum(count))
        w   = np.array([np.sum(p[:k]) for k in range(bins)])
        m   = np.array([np.dot(p[:k],bincenter[:k]) for k in range(bins)])
        mT  = np.dot(p,bincenter)
        
        sB  = np.array([(mT * w[k] - m[k])**2/(w[k]*(1.-w[k])) if w[k]*(1.-w[k]) > 0 else 0 for k in range(bins)])
        idx = np.argmax(sB)
        
        return bincenter[idx]




    def CellDivisionTrajectory(self,dataID, discretize_by = 'length', sisterdata = None):
        
        if not discretize_by in self.keylist_stripped:  raise KeyError('key not found')
        # not sure if this is needed
        if sisterdata is None:  sisterdata = self.__sisterdata
        if sisterdata:          keysuffix = ['A','B']    
        else:                   keysuffix = ['']
        
        # return two pandas-dataframes for the two trajectories usually contained in one file as two elements of a list
        # as the two sisters do not need to have the same number of cell divisions, a single dataframe might cause problems with variably lengthed trajectories
        ret = list()
        for ks in keysuffix:
            # use 'discretize_by' as the column name that should serve as measurement that indicates the division
            # in our data, both 'length' and 'cellsize' should be OK.
            # get list of changes between measurements. Usually this distribution is bimodal:
            #  * scattering around small values for normal data-points
            #  * scattering around 'large' negative values for cell divisions
            diffdata  = np.diff(self.__data[dataID][discretize_by + ks])
            # estimate a threshold from the data between the two peaks in the bimodal distribution with Otsu's method, then get indices of these transitions
            index_div = np.where(diffdata < self.otsu(diffdata))[0].flatten()
            # timepoint of division is assumed to be the average before and after the drop in signal
            time_div  = 0.5 * np.array(self.__data[dataID]['time' + ks][index_div + 1]) + 0.5 * np.array(self.__data[dataID]['time'+ks][index_div])
            
            # store results in dictionary, which can be easier transformed into pandas-dataframe
            # values computed for each cell cycle are
            #  * generation time
            #  * size at birth
            #  * size at division
            #  * Least Mean Squares Estimator for the (exponential) growth rate over the full division cycle
            ret_ks = dict()
            ret_ks['generationtime' + ks]          = np.diff(time_div)
            ret_ks[discretize_by + '_birth' + ks]  = np.array(self.__data[dataID][discretize_by + ks][index_div+1])[:-1]
            ret_ks[discretize_by + '_final' + ks]  = np.array(self.__data[dataID][discretize_by + ks][index_div])[1:]
            ret_ks['growth_' + discretize_by + ks] = np.array([
                            self.LMSQ(
                                self.__data[dataID][discretize_by + ks][index_div[i]+1:index_div[i+1]+1],           # x-values
                                np.log(self.__data[dataID][discretize_by + ks][index_div[i]+1:index_div[i+1]+1]),   # y-values
                                cov=False)[1] for i in range(len(index_div)-1)])
            
            # we have everything, now make a dataframe
            ret.append(pd.DataFrame(ret_ks))
        return ret

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




        
