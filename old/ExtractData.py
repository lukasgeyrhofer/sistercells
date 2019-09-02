#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import argparse
import sys,math
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default = [])
    parser.add_argument("-C","--readColumns",default=["timeA","lengthA"],nargs="*",type=str)
    parser.add_argument("-v","--verbose",default=False,action="store_true")
    args = parser.parse_args()
    
    for i,fn in enumerate(args.infiles):
        outfn = os.path.splitext(os.path.basename(fn))[0] + '.txt'
        data = pd.read_excel(fn)
        if args.verbose:
            print('{:d} \033[91m{:s}\033[0m'.format(i,fn))
        
        outdata = np.array([data[key] for key in args.readColumns if key in data.keys()]).T
        
        np.savetxt(outfn, outdata)

if __name__ == "__main__":
    main()
    
    
