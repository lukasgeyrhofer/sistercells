#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys,math
import argparse

import sistercellclass as scc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",     default = [], nargs = "*")
    parser.add_argument("-n", "--TrajCount",   default = 10, type = int)
    parser.add_argument("-o", "--baseoutname", default = "RandTraj", type = str)
    parser.add_argument("-F", "--ForceLength", default = False, action = "store_true")
    parser.add_argument("-v", "--verbose",     default = False, action = "store_true")
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    for i,rt in enumerate(data.RandomTrajectories(count = args.TrajCount, force_length = args.ForceLength)):
        fn = args.baseoutname + '{:03d}.xls'.format(i)
        if args.verbose: print("writing file '{}'".format(fn))
        rt.to_excel(fn, index = False)


if __name__ == "__main__":
    main()
    
    
