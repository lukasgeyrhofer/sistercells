#!/usr/bin/env python

import numpy as np
import pandas as pd

class SisterCellData(object):
    def __init__(self,**kwargs):
        self.__infiles = kwargs.get('infiles',[])
        
        self.__sisterdata = kwargs.get('sisterdata',True)
        
        self.__data = list()
        self.__dataorigin = list()
        self.__keylist = list()
        
        for filename in self.__infiles:
            tmpdata = pd.read_excel(filename)
            self.__data.append(tmpdata)
            self.__dataorigin.append(filename)
            for k in tmpdata.keys():
                if not str(k) in self.__keylist:
                    self.__keylist.append(str(k))
        
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
        if sisterdata is None:  sisterdata = self.__sisterdata
        if sisterdata:          keysuffix = ['A','B']    
        else:                   keysuffix = ['']
        
        ret = list()
        for ks in keysuffix:
            diffdata  = np.diff(self.__data[dataID][discretize_by + ks])
            index_div = np.where(diffdata < self.otsu(diffdata))[0].flatten()
            time_div  = 0.5 * np.array(self.__data[dataID]['time' + ks][index_div + 1]) + 0.5 * np.array(self.__data[dataID]['time'+ks][index_div])
            
            ret_ks = dict()
            ret_ks['generationtime' + ks]          = np.diff(time_div)
            ret_ks[discretize_by + '_birth' + ks]  = np.array(self.__data[dataID][discretize_by + ks][index_div+1])[:-1]
            ret_ks[discretize_by + '_final' + ks]  = np.array(self.__data[dataID][discretize_by + ks][index_div])[1:]
            ret_ks['growth_' + discretize_by + ks] = np.array([
                            self.LMSQ(
                                self.__data[dataID][discretize_by + ks][index_div[i]+1:index_div[i+1]+1],           # x-values
                                np.log(self.__data[dataID][discretize_by + ks][index_div[i]+1:index_div[i+1]+1]),   # y-values
                                cov=False)[1] for i in range(len(index_div)-1)])
            
            ret.append(pd.DataFrame(ret_ks))
        return ret


    def __getitem__(self,key):
        return self.__data[key]
    
    def __getattr__(self,key):
        if key == "filenames":
            return self.__dataorigin
        elif key == "keylist":
            return self.__keylist
        elif key == "keylist_stripped":
            return list(set([s.strip('AB ') for s in self.__keylist]))
    
    def __len__(self):
        return len(self.__data)

    def __iter__(self):
        dataIDs = np.arange(len(self),dtype=int)
        for dataID,origin,data in zip(dataIDs,self.__dataorigin,self.__data):
            yield dataID,origin,data




        
