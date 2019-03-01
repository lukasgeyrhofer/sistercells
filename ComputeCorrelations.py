#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse

import sistercellclass as ssc


def main():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i", "--infiles",       default = [],           nargs = "*")
    parser_io.add_argument("-o", "--OutFilePrefix", default = 'corrmatrix', type = str)
    parser_io.add_argument("-v", "--verbose",       default = False,        action = "store_true")
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-k","--DiscretizeKey", default = 'length',     type = str)
    parser_alg.add_argument("-m","--MaxLag",        default = 20,           type = int)
    parser_alg.add_argument("-S","--Symmetrize",    default = False,        action = "store_true")
    parser_alg.add_argument("-N","--Normalize",     default = False,        action = "store_true")
    args = parser.parse_args()


    data = ssc.SisterCellData(**vars(args))

    corrmatrix_sumAB   = dict()
    corrmatrix_sumA    = dict()
    corrmatrix_sumB    = dict()
    corrmatrix_count   = dict()
    
    
    for dataID,fn,x in data:
        if args.verbose: print(dataID,fn)
        
        trajA,trajB = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)

        for corrkey in [ k.strip('AB ') for k in trajA.keys() if k[:4].upper() != 'TIME']:
            if not corrkey in corrmatrix_sumAB.keys():
                corrmatrix_sumAB[corrkey] = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                corrmatrix_sumA[corrkey]  = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                corrmatrix_sumB[corrkey]  = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                corrmatrix_count[corrkey] = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
        
            for i in range(min(args.MaxLag,len(trajA))):
                for j in range(min(args.MaxLag,len(trajB))):
                    corrmatrix_sumAB[corrkey][i,j]   += trajA[corrkey + 'A'][i] * trajB[corrkey + 'B'][j]
                    corrmatrix_sumA[corrkey][i,j]    += trajA[corrkey + 'A'][i]
                    corrmatrix_sumB[corrkey][i,j]    += trajB[corrkey + 'B'][j]
                    corrmatrix_count[corrkey][i,j]   += 1
                

    # compute correlation & output
    cm = dict()
    for corrkey in corrmatrix_sumAB.keys():
        cm[corrkey] = (corrmatrix_sumAB[corrkey] - corrmatrix_sumA[corrkey] * corrmatrix_sumB[corrkey] / corrmatrix_count[corrkey])/corrmatrix_count[corrkey]

        fp = open(args.OutFilePrefix + '_' + corrkey,'w')
        for i in range(np.shape(cm[corrkey])[0]):
            for j in range(np.shape(cm[corrkey])[1]):
                outvalue = cm[corrkey][i,j]
                if args.Symmetrize: outvalue = 0.5*(cm[corrkey][i,j] + cm[corrkey][j,i])
                if args.Normalize:  outvalue /= cm[corrkey][0,0]
                fp.write('{:3d} {:3d} {:14.6e}\n'.format(i,j,outvalue))
            fp.write('\n')
        fp.close()
        


if __name__ == "__main__":
    main()

