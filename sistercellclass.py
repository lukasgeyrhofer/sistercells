#!/usr/bin/env python

import numpy as np
import pandas as pd

class SisterCellData(object):
    def __init__(self,**kwargs):
        self.__infiles = kwargs.get('infiles',[])
        
        self.__data = list()
        self.__dataorigin = list()
        
        for filename in self.__infiles:
            self.__data.append(pd.read_excel(filename))
            self.__dataorigin.append(filename)
        
        if not len(self) > 0:
            raise IOError("no data loaded")
        
    def __len__(self):
        return len(self.__data)

    def __iter__(self):
        for origin,data in zip(self.__dataorigin,self.__data):
            yield origin,data

