#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys,math
import numpy as np
import argparse


def LMSQ(x,y,cov = True):
    """
    Least Mean Squares estimator
    for a linear interpolation between x,y: y ~ a + b x
    Returns also covariance matrix for the estimated parameters a,b if 'cov' is True
    
    Algorithm:
    A = ( 1  ... 1  )
        ( x1 ... xn )
        
    y = ( y1 ... yn ).T
    
    p = (a b).T
    
    E[p]     = inv(A.T * A) * A.T * y
    Cov[p,p] = sigma2 * inv(A.T * A)
    sigma2   = E[ ( y - A*E[p] )^2 ]
    """
    
    n   = len(x)
    sx  = np.sum(x)
    sy  = np.sum(y)
    sxx = np.dot(x,x)
    sxy = np.dot(x,y)
    syy = np.dot(y,y)
    
    # estimate parameters
    denom     = (n*sxx-sx*sx)
    b         = (n*sxy - sx*sy)/denom
    a         = (sy-b*sx)/n
    estimate  = np.array([a,b],dtype=np.float)

    if cov:
        # estimate covariance matrix of estimated parameters
        sigma2    = syy + n*a*a + b*b*sxx + 2*a*b*sx - 2*a*sy - 2*b*sxy # variance of deviations from linear line
        covmatrix = sigma2 / denom * np.array([[sxx,-sx],[-sx,n]],dtype=np.float)

        return estimate,covmatrix
    else:
        return estimate


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


    diffcoeff , dc_cov = LMSQ(np.arange(mc), dt_var[:mc], cov = True)
    
    print("# {:14.6e} {:14.6e}".format(diffcoeff[1],np.sqrt(dc_cov[1,1])), file = sys.stderr)
    
    for i,v in enumerate(zip(dt_mean[:mc],dt_var[:mc],dt_count[:mc])):
        print('{:4d} {:14.6f} {:14.6f} {:.0f}'.format(i,*v))

if __name__ == "__main__":
    main()
