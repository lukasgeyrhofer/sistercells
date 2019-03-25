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
    args = parser.parse_args()


    data = ssc.SisterCellData(**vars(args))
    print(data)

    corrmatrix_sumAB   = dict()
    corrmatrix_sumA    = dict()
    corrmatrix_sumB    = dict()
    corrmatrix_count   = dict()
    
    if not args.ACFFilePrefix is None:
        acf_sumAA   = dict()
        acf_sumA    = dict()
        acf_countAA = dict()
        acf_countA  = dict()
    
    
    for dataID in range(len(data)):
        if args.verbose: print('{:4d} \033[91m{:s}\033[0m'.format(dataID,data.filenames[dataID]))
        
        cmAB,cmA,cmB,cmcount = data.lineagecorrelation(dataID, maxlen = args.MaxLag)
        
        for key in cmAB.keys():
            if not key in corrmatrix_sumAB.keys():
                corrmatrix_sumAB[key]  = cmAB[key]
                corrmatrix_sumA [key]  = cmA[key]
                corrmatrix_sumB [key]  = cmB[key]
                corrmatrix_count[key]  = cmcount[key]
            else:
                corrmatrix_sumAB[key] += cmAB[key]
                corrmatrix_sumA [key] += cmA[key]
                corrmatrix_sumB [key] += cmB[key]
                corrmatrix_count[key] += cmcount[key]

        if not args.ACFFilePrefix is None:
            
            acfAA, acfcAA, acfA, acfcA = data.autocorrelation_restricted(dataID, maxlen = args.MaxLag)
            
            for key in acfAA.keys():
                if not key in acf_sumAA.keys():
                    acf_sumAA  [corrkey] = np.zeros(args.MaxLag)
                    acf_countAA[corrkey] = np.zeros(args.MaxLag)
                    acf_sumA   [corrkey] = 0
                    acf_countA [corrkey] = 0
                    
                acf_sumAA  [corrkey][:len(acfAA[corrkey])] += acfAA[corrkey]
                acf_countAA[corrkey][:len(acfAA[corrkey])] += acfcAA[corrkey]
                acf_sumA   [corrkey] += acfA[corrkey]
                acf_countA [corrkey] += acfcA[corrkey]
                

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
                fp.write('{:3d} {:3d} {:14.6e}\n'.format(i+1,j+1,outvalue)) # add 1, since first data point is first (!) generation of sisters
            fp.write('\n')
        fp.close()
        
        if not args.ACFFilePrefix is None:
            acf[corrkey] = acf_sumAA[corrkey]/acf_countAA[corrkey] - acf_sumA[corrkey] * acf_sumA[corrkey] / ( acf_countA[corrkey] * acf_countA[corrkey] )
            if args.Normalize: acf[corrkey] = 1./acf[corrkey][0]
            np.savetxt(args.ACFFilePrefix + '_' + corrkey,np.array([np.arange(len(acf[corrkey])),acf[corrkey]]).T)

if __name__ == "__main__":
    main()

