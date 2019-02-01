#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import sistercellclass as scc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-k","--DiscretizeKey",default="length",type=str)
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    correlationfunctions = dict()
    
    for x in data:
        for y in data.CellDivisionTrajectory(x[0],discretize_by = args.DiscretizeKey):
            for k in y.keys():
                if k != 'time':
                    if not k in correlationfunctions.keys():
                        correlationfunctions[k] = list()
                    correlationfunctions[k].append(np.correlate(y[k],y[k]))
    cf_avg = dict()
    for k in correlationfunctions.keys():
        l = np.min([len(a) for a in correlationfunctions[k]])
        cf_avg[k] = np.avg([a[:l] for a in correlationfunctions[k]])
        print k
        print cf_avg
            
    
    
if __name__ == "__main__":
    main()
