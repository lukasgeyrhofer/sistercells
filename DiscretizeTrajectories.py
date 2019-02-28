#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

import sistercellclass as scc


def autocorrelation (x) :
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    """
    xp = x-np.mean(x)
    f = np.fft.fft(xp)
    p = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
    pi = np.fft.ifft(p)
    return np.real(pi)[:x.size/2]/np.sum(xp**2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*",default=[])
    parser.add_argument("-k","--DiscretizeKey",default="length",type=str)
    parser.add_argument("-o","--outfileprefix",default="ACF",type=str)
    parser.add_argument("-m","--minlength",default=10,type=int)
    parser.add_argument("-S","--averageSisters",default=False,action="store_true")
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    correlationfunctions = dict()
    
    for x in data:
        for y in data.CellDivisionTrajectory(x[0],discretize_by = args.DiscretizeKey):
            for k in y.keys():
                if k != 'time' and len(y[k]) > args.minlength:
                    k0 = k
                    if args.averageSisters:
                        k0 = k.rstrip('AB')
                    if not k0 in correlationfunctions.keys():
                        correlationfunctions[k0] = list()
                    correlationfunctions[k0].append(autocorrelation(y[k]))
                    
    for k in correlationfunctions.keys():
        maxL = np.max([len(a) for a in correlationfunctions[k]])
        cf_sum   = np.zeros(maxL,dtype = np.float)
        cf_count = np.zeros(maxL,dtype = np.float)
        for acf in correlationfunctions[k]:
            if not np.any(np.isnan(acf)):
                cf_sum[:len(acf)] += acf
                cf_count[:len(acf)] += 1
        if cf_count[-1] == 0:
            maxL = np.argmin(cf_count > 0)
        np.savetxt(args.outfileprefix + "_" + k.replace(" ","-"),np.array([np.arange(maxL),cf_sum[:maxL]/cf_count[:maxL]]).T)
            
    
    
if __name__ == "__main__":
    main()
