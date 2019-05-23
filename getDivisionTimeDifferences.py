#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse
import sys,math
import os

import sistercellclass as ssc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs= "*", default = [])
    parser.add_argument("-o", "--outfilesuffix", default = "dtd", type = str)
    parser.add_argument("-v", "--verbose", default = False, action = "store_true")
    parser.add_argument("-O", "--Get_OTSU_From_Both", default = True, action = "store_false")
    args = parser.parse_args()

    data = ssc.SisterCellData(**vars(args))
    
    for i in range(len(data)):
        if args.verbose: print(data.filenames[i])
        print(i)
        dtd = data.DivisionTimeDifferences(i)
        outfilename = os.path.splitext(os.path.basename(data.filenames[i]))[0] + '.' + args.outfilesuffix
        np.savetxt(outfilename, np.array([np.arange(len(dtd)),dtd]).T)


if __name__ == "__main__":
    main()
