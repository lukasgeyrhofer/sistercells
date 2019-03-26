#!/usr/bin/env python3

from __future__ import print_function

import numpy as np
import argparse
import sys,math


#def autocorrelation(x, normalize = False) :
    #"""
    #Compute the autocorrelation of the signal, based on the properties of the
    #power spectral density of the signal.
    #"""
    #xp  = x-np.mean(x)
    #f   = np.fft.fft(xp)
    #p   = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
    #pi  = np.fft.ifft(p)
    #acf = np.real(pi)[:x.size/2]
    #if normalize: acf *= 1./np.sum(xp**2)
    #return acf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",     default = [], type=str, nargs = "*")
    parser.add_argument("-o", "--outfile",     default = None)
    parser.add_argument("-T", "--timecolumn",  default = 0, type = int)
    parser.add_argument("-D", "--datacolumn",  default = 1, type = int)
    parser.add_argument("-M", "--maxtime",     default = None, type = int)
    parser.add_argument("-N", "--normalize",   default = False, action = "store_true")
    parser.add_argument("-t", "--transitive",  default = True, action = "store_false")
    args = parser.parse_args()

    
    dt          = list()
    acf_sumAA   = list()
    acf_countAA = np.array([],dtype=np.float)
    acf_sumA    = 0
    acf_countA  = 0
    
    for fn in args.infiles:
        try:
            data = np.genfromtxt(fn)
        except:
            continue
        
        maxtime = np.shape(data)[0]
        
        t = data[:,args.timecolumn]
        x = data[:,args.datacolumn]
        
        dt.append(np.mean(np.diff(t)))
        
        if len(acf_sumAA) < maxtime:
            acf_sumAA   = np.concatenate([acf_sumAA,  np.zeros(maxtime - len(acf_sumAA)   + 1)])
            acf_countAA = np.concatenate([acf_countAA,np.zeros(maxtime - len(acf_countAA) + 1)])

        for i in range(maxtime):
            if args.transitive:
                for j in range(maxtime - i):
                    acf_sumAA[i]       += x[j] * x[j+i]
                    acf_countAA[i] += 1
            else:
                acf_sumAA[i]       += x[0] * x[i]
                acf_countAA[i] += 1
        
            acf_sumA   += x[i]
            acf_countA += 1
        
    for cdt in dt:
        if not math.isclose(cdt,dt[0]):
            raise ValueError('Data files have different time steps')

    timedigits = int(np.ceil(-np.log10(dt[0])))

    maxtime = len(acf_sumAA)
    if not args.maxtime is None:
        if maxtime > args.maxtime:
            maxtime = args.maxtime

    acf = acf_sumAA[:maxtime] / acf_countAA[:maxtime] - (acf_sumA/acf_countA)**2
    
    if args.normalize:
        acf /= acf[0]
    
    try:    fp = open(args.outfile,'w')
    except: fp = sys.stdout
        
    for values in zip(np.arange(maxtime) * dt[0],acf,acf_countAA):
        fp.write("{:.{}f} {:14.6e} {}\n".format(values[0],timedigits,values[1],values[2]))
    
    fp.close()
    
    
if __name__ == "__main__":
    main()
