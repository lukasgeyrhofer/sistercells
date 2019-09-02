#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import pandas as pd

import sistercellclass as scc

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-k","--key",default="length",type=str)
    parser.add_argument("-o","--outfile",default="out",type=str)
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    if not args.key in data.keylist_stripped:
        raise KeyError
    
    fig,axes = plt.subplots(figsize = (20,10), dpi = 100, ncols = 2, nrows = len(data), sharey = True, sharex = True)
    
    for dID,fn,x in data:
        dA = np.diff(x[args.key + 'A'])
        dB = np.diff(x[args.key + 'B'])
        axes[dID,0].hist(dA,log=True,bins=50)
        axes[dID,0].axvline(x = data.otsu(dA),c='red')
        axes[dID,1].hist(dB,log=True,bins=50)
        axes[dID,1].axvline(x = data.otsu(dB),c='red')
        
    fig.savefig(args.outfile + '.svg')


if __name__ == "__main__":
    main()
