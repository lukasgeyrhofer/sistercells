#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse
import sys,math
import os
import sistercellclass as scc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-o","--outfilesuffix",default="discretized")
    parser.add_argument("-D", "--DiscretizeKey", default = "length", type   = str)
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    for dataID,fn,x in data:
        trajA,trajB = data.CellDivisionTrajectory(dataID, discretize_by = args.DiscretizeKey)
        trajA.to_csv(os.path.basename(fn)[:-3] + args.outfilesuffix + 'A',sep = ' ',index_label='# generation')
        trajB.to_csv(os.path.basename(fn)[:-3] + args.outfilesuffix + 'B',sep = ' ',index_label='# generation')

if __name__ == "__main__":
    main()
