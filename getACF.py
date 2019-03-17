#!/usr/bin/env python

from __future__ import print_function

import argparse
import numpy as np
import sys,math
import pandas as pd

import sistercellclass as scc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",       default = [],       nargs  = "*")
    parser.add_argument("-o", "--outfileprefix", default = 'ACF',    type   = str)
    parser.add_argument("-M", "--maxtime",       default = None,     type   = int)
    parser.add_argument("-D", "--DiscretizeKey", default = "length", type   = str)
    parser.add_argument("-v", "--verbose",       default = False,    action = "store_true")
    parser.add_argument("-N", "--normalize",     default = False,    action = "store_true")
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    acf_AB      = dict()
    acf_ABcount = dict()
    acf_A       = dict()
    acf_Acount  = dict()
    
    for dataID,fn,x in data:
        if args.verbose: print(fn)
        traj = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)
        
        for a in range(2):
            for k in data.keylist_stripped:
                trajlen = np.max([args.maxtime,len(traj[a][k])])
                if not k in acf_AB.keys():
                    acf_AB[k]      = np.array([],dtype=np.float)
                    acf_ABcount[k] = np.array([],dtype=np.float)
                    
                    acf_A[k]       = np.array([],dtype=np.float)
                    acf_Acount[k]  = np.array([],dtype=np.float)

                for i in range(trajlen):
                    if len(acf_A[k]) > i:
                        acf_A[k]      = np.concatenate([acf_A[k],      np.zeros(1)])
                        acf_Acount[k] = np.concatenate([acf_Acount[k], np.zeros(1)])
                    acf_A[k][i]      += traj[a][k][i]
                    acf_Acount[k][i] += 1.
                    
                    for j in range(trajlen - i):
                        if len(acf_AB[k]) > j:
                            acf_AB[k]      = np.concatenate([acf_AB[k],      np.zeros(1)])
                            acf_ABcount[k] = np.concatenate([acf_ABcount[k], np.zeros(1)])
                        acf_AB[k][j]      += traj[a][k][i] * traj[a][k][i+j]
                        acf_ABcount[k][j] += 1.
                
        
    for k in acf_AB.keys():
        assert len(acf_AB) == len(acf_A)
        
        tmpacf = acf_AB[k]/acf_ABcount[k] - (acf_A[k]/acf_Acount[k])**2
        
        if args.normalize:
            tmpacf /= tmpacf[0]
            
        t = np.arange(len(tmpacf))
        np.savetxt(args.outfileprefix + '_' + k.replace(' ','-'),np.array([t,tmpacf]).T)



if __name__ == "__main__":
    main()

    
