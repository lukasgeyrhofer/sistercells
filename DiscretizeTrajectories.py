#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import sistercellclass as scc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-k","--DiscretizeKey",default="lenght",type=str)
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    for x in data:
        for y in data.CellDivisionTrajectory(x[0],discretize_by = args.DiscretizeKey):
            print y
            print
            
    
    
if __name__ == "__main__":
    main()
