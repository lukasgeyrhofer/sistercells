#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    args = parser.parse_args()

    for filename in args.infiles:
        tmpdata = pd.read_excel(filename)
        keys = np.array(tmpdata.keys(),dtype=str)
        tracelength = len(tmpdata['timeA'])
        print '\033[91m{:40s}\033[0m \033[92m{:4d}\033[0m  ['.format(filename,tracelength) + ', '.join(keys) + ']'

if __name__ == "__main__":
    main()
