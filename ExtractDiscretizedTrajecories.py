#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse
import sys,math
import os
import sistercellclass as scc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs = "*", default = [])
    parser.add_argument("-o", "--outfilesuffix", default = "discretized", type = str)
    parser.add_argument("-D", "--DiscretizeKey", default = "length", type = str)
    parser.add_argument("-H", "--Histograms", default = False, action = "store_true")
    parser.add_argument("-b", "--bins", default=100,type=int)
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    if args.Histograms:
        histodata = dict()
    
    for dataID,fn,x in data:
        trajA,trajB = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)
        trajA.to_csv(os.path.basename(fn)[:-3] + args.outfilesuffix + 'A',sep = ' ',index_label='# generation')
        trajB.to_csv(os.path.basename(fn)[:-3] + args.outfilesuffix + 'B',sep = ' ',index_label='# generation')
        
        if args.Histograms:
            for k in trajA.keys():
                if not k in histodata.keys():
                    histodata[k] = trajA[k]
                else:
                    histodata[k] = np.concatenate([histodata[k],trajA[k]])
                histodata[k] = np.concatenate([histodata[k],trajB[k]])
                

    if args.Histograms:
        for k in histodata.keys():
            h,b = np.histogram(histodata[k],bins = args.bins)
            b = b[:-1] + 0.5 * np.diff(b)
            
            np.savetxt(k.replace(' ','-') + '.histo',np.array([b,h]).T)
                                               
        

if __name__ == "__main__":
    main()
