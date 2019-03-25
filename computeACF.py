#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import argparse
import sys


def autocorrelation(x, normalize = False) :
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    """
    xp  = x-np.mean(x)
    f   = np.fft.fft(xp)
    p   = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
    pi  = np.fft.ifft(p)
    acf = np.real(pi)[:x.size/2]
    if normalize: acf *= 1./np.sum(xp**2)
    return acf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile",type=str)
    parser.add_argument("-o","--outfile",default=None)
    parser.add_argument("-T","--timecolumn",default=0,type=int)
    parser.add_argument("-D","--datacolumn",default=1,type=int)
    parser.add_argument("-M","--maxtime",default=None,type=int)
    parser.add_argument("-N","--normalize",default=False,action="store_true")
    args = parser.parse_args()
    
    try:    data = np.genfromtxt(args.infile)
    except: raise IOError("could not open file '{}'".format(args.infile))
    
    maxtime = np.shape(data)[0]
    if not args.maxtime is None:
        maxtime = np.min(maxtime,args.maxtime)
    
    t = data[:maxtime,args.timecolumn]
    x = data[:maxtime,args.datacolumn]
    a = autocorrelation(x, normalize = args.normalize)
    
    try:    fp = open(args.outfile,'w')
    except: fp = sys.stdout
    
    dt = np.mean(np.diff(t))
    timedigits = int(np.ceil(-np.log10(dt)))
    
    for values in zip(np.arange(maxtime) * dt,a):
        fp.write("{:.{}f} {:14.6e}\n".format(values[0],timedigits,values[1]))
    
    fp.close()
    
    
if __name__ == "__main__":
    main()
