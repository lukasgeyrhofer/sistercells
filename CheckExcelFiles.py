#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    args = parser.parse_args()

    maxinfilelength = np.max([len(fn) for fn in args.infiles])

    for filename in args.infiles:
        tmpdata = pd.read_excel(filename)
        keys = np.array(tmpdata.keys(),dtype=str)
        tracelength = len(tmpdata[keys[0]])
        print '\033[91m{:{}s}\033[0m \033[92m{:4d}\033[0m  ['.format(filename,maxinfilelength,tracelength) + ', '.join(keys) + ']'

if __name__ == "__main__":
    main()
