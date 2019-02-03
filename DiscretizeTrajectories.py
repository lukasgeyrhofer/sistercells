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
    args = parser.parse_args()
    
    data = scc.SisterCellData(**vars(args))
    
    correlationfunctions = dict()
    
    for x in data:
        for y in data.CellDivisionTrajectory(x[0],discretize_by = args.DiscretizeKey):
            for k in y.keys():
                if k != 'time':
                    if not k in correlationfunctions.keys():
                        correlationfunctions[k] = list()
                    correlationfunctions[k].append(np.correlate(y[k],y[k]))
    cf_avg = dict()
    for k in correlationfunctions.keys():
        l = np.min([len(a) for a in correlationfunctions[k]])
        cf_avg[k] = np.avg([a[:l] for a in correlationfunctions[k]])
        print k
        print cf_avg
            
    
    
if __name__ == "__main__":
    main()
