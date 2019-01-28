#!/usr/bin/env python


import numpy as np
import pandas as pd
import argparse
import sys,math

import sistercellclass as scc

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    args=parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    fig,axes  = plt.subplots(nrows=len(data),ncols=1,figsize=(10, 100), dpi=100, frameon = False, sharex = True, sharey = True)
    
    maxtimes = list()
    
    
    for i,d in enumerate(data):
        maxtimes.append(np.max(d[1][u'timeA']))
        print 'open file \033[91m{:70s}\033[0m maxtime: {:5.2f}'.format(d[0],np.max(d[1][u'timeA']))
        axes[i].plot(d[1][u'timeA'],d[1][u'lengthA'])
        axes[i].plot(d[1][u'timeB'],d[1][u'lengthB'])
    
    fig.savefig('alltrajectories.svg')
    
    fig2 = plt.figure()
    plt.hist(maxtimes,range=(0,100),bins=101)
    fig2.savefig('maxtimehist.svg')

if __name__ == "__main__":
    main()
    
