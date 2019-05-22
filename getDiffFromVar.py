#!/usr/bin/env python3

import sys,math
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",  default = [],   nargs = "*",)
    parser.add_argument("-m", "--MinCount", default = None, type = int)
    args = parser.parse_args()

    dt_sum   = np.zeros(1, dtype = np.float)
    dt2_sum  = np.zeros(1, dtype = np.float)
    dt_count = np.zeros(1, dtype = np.float)

    for fn in args.infiles:
        try:
            data = np.genfromtxt(fn)
        except:
            continue
        
        n = len(data)
        if n > len(dt_sum):
            dt_sum   = np.concatenate([dt_sum,   np.zeros(n - len(dt_sum),   dtype = np.float)])
            dt2_sum  = np.concatenate([dt2_sum,  np.zeros(n - len(dt2_sum),  dtype = np.float)])
            dt_count = np.concatenate([dt_count, np.zeros(n - len(dt_count), dtype = np.float)])
        
        dt_sum[:n]   += data[:,1]
        dt2_sum[:n]  += data[:,1] * data[:,1]
        dt_count[:n] += np.ones(n)
    
    dt_mean = dt_sum / dt_count
    dt_var  = dt2_sum / dt_count - dt_mean ** 2
    
    if not args.MinCount is None:
        mc = min(np.where(dt_count < args.MinCount)[0])
    else:
        mc = len(dt_mean)
        
    for i,v in enumerate(zip(dt_mean[:mc],dt_var[:mc],dt_count[:mc])):
        print('{:4d} {:14.6f} {:14.6f} {:.0f}'.format(i,*v))

if __name__ == "__main__":
    main()
