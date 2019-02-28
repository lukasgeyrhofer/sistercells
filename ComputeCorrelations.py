#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse
import sys,math

import sistercellclass as ssc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",default=[],nargs="*")
    parser.add_argument("-k","--DiscretizeKey",default="length",type=str)
    parser.add_argument("-K","--CorrelationKey",default="length_final",type=str)
    parser.add_argument("-m","--MaxLag",type=int,default=20)
    parser.add_argument("-o","--outfile",type=str,default='corrmatrix.txt')
    parser.add_argument("-v","--verbose",default=False,action="store_true")
    parser.add_argument("-S","--symmetrize",default=False,action="store_true")
    args = parser.parse_args()


    data = ssc.SisterCellData(**vars(args))

    corrmatrix_sumAB   = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
    corrmatrix_sumA    = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
    corrmatrix_sumB    = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
    
    corrmatrix_count   = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
    
    avgA_sum         = 0.
    avgB_sum         = 0.
    avgA_count       = 0.
    avgB_count       = 0.


    for dataID,fn,x in data:
        if args.verbose:
            print(dataID,fn)
        
        trajA,trajB = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)
        
        for i in range(min(args.MaxLag,len(trajA))):
            for j in range(min(args.MaxLag,len(trajB))):
                corrmatrix_sumAB[i,j]   += trajA[args.CorrelationKey + 'A'][i] * trajB[args.CorrelationKey + 'B'][j]
                corrmatrix_sumA[i,j]    += trajA[args.CorrelationKey + 'A'][i]
                corrmatrix_sumB[i,j]    += trajB[args.CorrelationKey + 'B'][j]
                
                corrmatrix_count[i,j]   += 1
                
    # compute correlation
    cm = (corrmatrix_sumAB - corrmatrix_sumA * corrmatrix_sumB / corrmatrix_count)/corrmatrix_count

    # symmetrize
    if args.symmetrize: cm = 0.5*(cm + cm.T)

    # output
    fp = open(args.outfile,'w')
    for i in range(np.shape(cm)[0]):
        for j in range(np.shape(cm)[1]):
            fp.write('{:3d} {:3d} {:14.6e}\n'.format(i,j,cm[i,j]))
        fp.write('\n')
    fp.close()
        


if __name__ == "__main__":
    main()

