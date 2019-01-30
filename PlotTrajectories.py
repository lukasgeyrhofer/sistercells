#!/usr/bin/env python


import numpy as np
import pandas as pd
import argparse
import sys,math

import sistercellclass as scc

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",         default = [], nargs = "*")
    parser.add_argument("-o","--outfile",         default = "alltrajectories", type = str)
    parser.add_argument("-k","--key",             default = "length", type = str)
    parser.add_argument("-H","--HistoTrajLength", default = False, action = "store_true")
    parser.add_argument("-v","--verbose",         default = False, action = "store_true")
    args=parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    assert args.key in data.keylist_stripped, 'option -k "KEY" only allows values in [' + ', '.join(data.keylist_stripped) + ']'
    
    fig,axes  = plt.subplots(nrows=len(data),ncols=1,figsize=(10, 100), dpi=100, frameon = False, sharex = True, sharey = True)
    
    maxtimes = list()
    
    for i,fn,d in data:
        maxtimes.append(np.max(d[u'timeA']))
        if args.verbose:
            print 'open file \033[91m{:s}\033[0m maxtime: {:5.2f}'.format(fn,np.max(d[u'timeA']))
        axes[i].plot(d[u'timeA'],d[args.key + 'A'])
        axes[i].plot(d[u'timeB'],d[args.key + 'B'])
    
    if args.verbose:
        print 'writing trajectories to \033[94m{:s}\033[0m'.format(args.outfile + '.svg')
    fig.savefig(args.outfile + '.svg')
    
    if args.HistoTrajLength:
        if args.verbose:
            print 'writing histogram to \033[94m{:s}\033[0m'.format('maxtimehist.svg')
        fig2 = plt.figure()
        plt.hist(maxtimes,range=(0,100),bins=101)
        fig2.savefig('maxtimehist.svg')

if __name__ == "__main__":
    main()
    
