#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys,math

import sistercellclass as scc

import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt


def otsu(x,historange = None, bins = None):
    if historange is None:
        if bins is None:
            count,binedges = np.histogram(x)
        else:
            count,binedges = np.histogram(x,bins=bins)
    else:
        if bins is None:
            count,binedges = np.histogram(x,range=historange)
        else:
            count,binedges = np.histogram(x,range=historange,bins=bins)
    bincenter = binedges[:-1] + .5 * np.diff(binedges)

    p   = count/float(sum(count))
    w   = np.array([np.sum(p[:k]) for k in range(bins)])
    m   = np.array([np.dot(p[:k],bincenter[:k]) for k in range(bins)])
    mT  = np.dot(p,bincenter)
    
    sB  = np.array([(mT * w[k] - m[k])**2/(w[k]*(1.-w[k])) if w[k]*(1.-w[k]) > 0 else 0 for k in range(bins)])
    idx = np.argmax(sB)
    
    return bincenter[idx]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-o","--outfile",default="diffhisto",type=str)
    args=parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))

    fig,axes  = plt.subplots(nrows=len(data),ncols=4,figsize=(10, 10), dpi=100, frameon = False, sharex = True, sharey = True)
    plt.yscale('log')
    
    for i,fn,d in data:
        bins = max(len(d['lengthA']) // 50,25)
        dlA = np.diff(d['lengthA'])
        dlB = np.diff(d['lengthB'])
        daA = np.diff(d['cellareaA'])
        daB = np.diff(d['cellareaB'])
        axes[i,0].hist(dlA,bins=bins)
        axes[i,0].axvline(x = otsu(dlA,bins=bins),c='red')
        axes[i,1].hist(dlB,bins=bins)
        axes[i,1].axvline(x = otsu(dlB,bins=bins),c='red')
        axes[i,2].hist(daA,bins=bins)
        axes[i,2].axvline(x = otsu(daA,bins=bins),c='red')
        axes[i,3].hist(daB,bins=bins)
        axes[i,3].axvline(x = otsu(daB,bins=bins),c='red')
        
        
        idx_lA = np.where(dlA < otsu(dlA,bins=bins))[0].flatten()
        idx_lB = np.where(dlB < otsu(dlB,bins=bins))[0].flatten()
        idx_aA = np.where(daA < otsu(daA,bins=bins))[0].flatten()
        idx_aB = np.where(daB < otsu(daB,bins=bins))[0].flatten()
        
        for t,l in zip(d['timeA'],d['lengthA']):
            print >> sys.stderr,t,l
        
        for i in idx_lA:
            print d['timeA'][i],d['lengthA'][i]
            #print d['timeA'][i+1],d['lengthA'][i+1]
        exit(1)
        
    fig.savefig(args.outfile + '.svg')
    
    
    
if __name__ == "__main__":
    main()
