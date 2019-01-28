#!/usr/bin/env python

import numpy as np
import pandas as pd

class SisterCellData(object):
    def __init__(self,**kwargs):
        self.__infiles = kwargs.get('infiles',[])
        
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
            raise IOError("no data loaded")

    def CellDivisionTrajectory(self,dataID, discretize_by = 'length', sisters = True):
        return None





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




        
