#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import argparse
import sys,math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    args = parser.parse_args()

    sumD   = 0
    sumDD  = 0
    countD = 0
    
    for fn in args.infiles:
        try:    data = np.genfromtxt(fn)
        except: continue

        g       = data[:,0] # generation
        dtd     = data[:,1] # division time difference
        
        sumD   += np.sum(dtd[1:] * dtd[1:] / g[1:])
        sumDD  += np.sum(dtd[1:] * dtd[1:] * dtd[1:] * dtd[1:] / (g[1:] * g[1:]))
        countD += len(g) - 1
    
    
    avgD = sumD/countD
    varD = sumDD/countD - avgD**2
    print('{:.6e} {:.6e} {}'.format(avgD,varD,countD))


if __name__ == "__main__":
    main()
