#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse
import sys,math

import sistercellclass as ssc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",default=[],nargs="*")
    parser.add_argument("-K","--DiscretizeKey",default='length',type=str)
    parser.add_argument("-N","--normalize",default=False,action="store_true")
    args = parser.parse_args()
    
    data = ssc.SisterCellData(**vars(args))
    
    histodata = np.array([],dtype=int)
    dt = data.timestep
    idt = 1./dt
    histo = np.zeros(0,dtype=int)
    
    for dataID,fn,x in data:
        trajA,trajB = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)
        
        
        for intdivtime in np.array(trajA['generationtime'] / dt,dtype=int):
            if intdivtime >= len(histo):
                histo = np.concatenate([histo,np.zeros(intdivtime - len(histo) + 1,dtype=int)])
            histo[intdivtime] += 1
    
        for intdivtime in np.array(trajB['generationtime'] / dt,dtype=int):
            if intdivtime >= len(histo):
                histo = np.concatenate([histo,np.zeros(intdivtime - len(histo) + 1,dtype=int)])
            histo[intdivtime] += 1
        
        
    if args.normalize:
        isum = 1./np.sum(histo)
        for t,h in enumerate(histo):
            print('{:.3f} {:.6e}'.format(t*dt,isum*h))
    else:
        for t,h in enumerate(histo):
            print('{:.3f} {:4d}'.format(t*dt,h))
        



if __name__ == "__main__":
    main()
    
    
    
