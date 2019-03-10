#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse

import sistercellclass as ssc


def main():
    parser = argparse.ArgumentParser()
    parser_io = parser.add_argument_group(description = "==== I/O parameters ====")
    parser_io.add_argument("-i", "--infiles",             default = [],           nargs = "*")
    parser_io.add_argument("-o", "--OutFilePrefix",       default = 'corrmatrix', type = str)
    parser_io.add_argument("-A", "--ACFFilePrefix",       default = None,         type = str)
    parser_io.add_argument("-v", "--verbose",             default = False,        action = "store_true")
    
    parser_alg = parser.add_argument_group(description = "==== Algorithm parameters ====")
    parser_alg.add_argument("-k","--DiscretizeKey",       default = 'length',     type = str)
    parser_alg.add_argument("-m","--MaxLag",              default = 10,           type = int)
    parser_alg.add_argument("-S","--Symmetrize",          default = False,        action = "store_true")
    parser_alg.add_argument("-N","--Normalize",           default = False,        action = "store_true")
    parser_alg.add_argument("-V","--NormVariance",        default = False,        action = "store_true")
    args = parser.parse_args()


    data = ssc.SisterCellData(**vars(args))

    corrmatrix_sumAB   = dict()
    corrmatrix_sumA    = dict()
    corrmatrix_sumB    = dict()
    corrmatrix_count   = dict()
    
    overall_sumAB = 0.
    overall_sumA  = 0.
    overall_sumB  = 0.
    overall_count = 0.
    
    if not args.ACFFilePrefix is None:
        acf_sumAB = dict()
        acf_sumA  = dict()
        acf_count = dict()
    
    for dataID,fn,x in data:
        if args.verbose: print('{:4d} \033[91m{:s}\033[0m'.format(dataID,fn))
        
        trajA,trajB = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)

        for corrkey in trajA.keys():
            if not corrkey in corrmatrix_sumAB.keys():
                corrmatrix_sumAB[corrkey] = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                corrmatrix_sumA [corrkey] = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                corrmatrix_sumB [corrkey] = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                corrmatrix_count[corrkey] = np.zeros((args.MaxLag,args.MaxLag),dtype=np.float)
                
                acf_sumAB[corrkey] = np.zeros(args.MaxLag + 1,dtype=np.float)
                acf_sumA [corrkey] = np.zeros(args.MaxLag + 1,dtype=np.float)
                acf_count[corrkey] = np.zeros(args.MaxLag + 1,dtype=np.float)
        
            for i in range(min(args.MaxLag,len(trajA))):
                for j in range(min(args.MaxLag,len(trajB))):
                    corrmatrix_sumAB[corrkey][i,j] += trajA[corrkey][i] * trajB[corrkey][j]
                    corrmatrix_sumA [corrkey][i,j] += trajA[corrkey][i]
                    corrmatrix_sumB [corrkey][i,j] += trajB[corrkey][j]
                    corrmatrix_count[corrkey][i,j] += 1
                    
                    overall_sumAB += trajA[corrkey][i] * trajB[corrkey][j]
                    overall_sumA  += trajA[corrkey][i]
                    overall_sumB  += trajB[corrkey][j]
                    overall_count += 1.
            
            if not args.ACFFilePrefix is None:
                for i in range(len(trajA)):
                    for j in range(np.min([args.MaxLag + 1,len(trajA) - i])):
                        acf_sumAB[corrkey][j] += trajA[corrkey][i] * trajA[corrkey][i+j]
                        acf_sumA[corrkey][j]  += trajA[corrkey][i]
                        acf_count[corrkey][j] += 1
                for i in range(len(trajB)):
                    for j in range(np.min([args.MaxLag + 1,len(trajB) - i])):
                        acf_sumAB[corrkey][j] += trajB[corrkey][i] * trajB[corrkey][i+j]
                        acf_sumA[corrkey][j]  += trajB[corrkey][i]
                        acf_count[corrkey][j] += 1
                        
                            
        
        
                

    if args.NormVariance:
        var = (overall_sumAB - overall_sumA * overall_sumB / overall_count) / overall_count

    # compute correlation & output
    cm = dict()
    acf = dict()
    for corrkey in corrmatrix_sumAB.keys():
        cm[corrkey] = (corrmatrix_sumAB[corrkey] - corrmatrix_sumA[corrkey] * corrmatrix_sumB[corrkey] / corrmatrix_count[corrkey])/corrmatrix_count[corrkey]

        fp = open(args.OutFilePrefix + '_' + corrkey,'w')
        for i in range(np.shape(cm[corrkey])[0]):
            for j in range(np.shape(cm[corrkey])[1]):
                outvalue = cm[corrkey][i,j]
                if args.Symmetrize:     outvalue  = 0.5*(cm[corrkey][i,j] + cm[corrkey][j,i])
                if args.Normalize:      outvalue /= cm[corrkey][0,0]
                elif args.NormVariance: outvalue /= var
                fp.write('{:3d} {:3d} {:14.6e}\n'.format(i+1,j+1,outvalue)) # add 1, since first data point is first (!) generation of sisters
            fp.write('\n')
        fp.close()
        
        if not args.ACFFilePrefix is None:
            acf[corrkey] = (acf_sumAB[corrkey] - acf_sumA[corrkey] * acf_sumA[corrkey] / acf_count[corrkey])/acf_count[corrkey]
            norm = 1.
            if args.Normalize: norm = 1./acf[corrkey][0]
            acf[corrkey] *= norm
            fp = open(args.ACFFilePrefix + '_' + corrkey,'w')
            for i in range(len(acf_sumAB[corrkey])):
                fp.write('{:3d} {:14.6e}\n'.format(i,acf[corrkey][i]))
            fp.close()

if __name__ == "__main__":
    main()

