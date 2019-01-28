#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
import sys,math

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    args=parser.parse_args()
    
    data = SisterCellData(**vars(args))
    
    plt.axis(xmin=0,xmax=10,ymin=0,ymax=5)
    fig  = plt.figure(figsize=(10, 100), dpi=100)
    grid = gridspec.GridSpec(ncols = 1, nrows = len(data), figure = fig)
    plots = list()
    
    maxtimes = list()
    
    
    for i,d in enumerate(data):
        maxtimes.append(np.max(d[1][u'timeA']))
        print 'open file \033[91m{:70s}\033[0m maxtime: {:5.2f}'.format(d[0],np.max(d[1][u'timeA']))
        plots.append(fig.add_subplot(grid[i,:]))
        plots[i].plot(d[1][u'timeA'],d[1][u'cellareaA'])
        plots[i].plot(d[1][u'timeB'],d[1][u'cellareaB'])
    
    fig.savefig('tmp.svg',)
    
    fig2 = plt.figure()
    plt.hist(maxtimes,range=(0,100),bins=101)
    fig2.savefig('maxtimehist.svg')

if __name__ == "__main__":
    main()
