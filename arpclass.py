#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math
import pandas as pd


class ARP(object):
    """
    2-dimensional autoregressive process with 2 noise terms
    """
    def __init__(self,**kwargs):
        # construct matrix A from eigenvalues and angle (0...1) between eigenvectors, assume first direction is (0,1)
        self.__A_eigval      = np.array([self.interval(x) for x in kwargs.get('eigval',[0.9,0.5])],dtype=np.float)
        self.__A_angle       = self.interval(kwargs.get('A_angle',0.33)) # angles between 0 .. 1
        evec1                = np.array([0,1],dtype=np.float)
        evec2                = np.dot(self.rotation(self.__A_angle),evec1)
        self.__A_eigvec      = (evec1 / np.linalg.norm(evec1), evec2 / np.linalg.norm(evec2))
        evmat                = np.array(self.__A_eigvec,dtype=float).T
        self.A               = np.matmul(np.matmul(evmat,np.diag(self.__A_eigval)), np.linalg.inv(evmat))
        
        # projection vector is also with respect to first EVec
        self.__alpha_angle   = self.interval(kwargs.get('alpha_angle',0.67))
        self.alpha           = np.dot(self.rotation(self.__alpha_angle),evec1)

        # noise, env
        self.noiseamplitudes = np.array(kwargs.get('noiseamplitudes',[1/np.sqrt(2.),1/np.sqrt(2.)]),dtype=np.float)

        self.experimenttype  = kwargs.get('experimenttype','sisters')
        self.default_steps   = kwargs.get('steps',10)
        
        # initialize all dynamical variables to start
        self.reset()
        
        
        # initialize analytical computations
        self.sumInf_AmATm   = self.compute_sumInf_AmATm()
        self.vardiff        = np.array([0],dtype=np.float)


    def random(self,mean = 0, sqrt_var = 1):
        return np.random.normal(loc = mean, scale = sqrt_var, size=2)


    def step(self):
        # noiseamplitudes = noise 0, env 1
        noiseA          = self.noiseamplitudes[0] * self.random()
        noiseB          = self.noiseamplitudes[0] * self.random()
        if self.experimenttype == 'sisters' or self.experimenttype == 'nonsisters':
            xiE         = self.random()
            noiseA     += self.noiseamplitudes[1] * xiE
            noiseB     += self.noiseamplitudes[1] * xiE
        else:
            noiseA     += self.noiseamplitudes[1] * self.random()
            noiseB     += self.noiseamplitudes[1] * self.random()

        xA_new = np.dot(self.A,self.xA[-1]) + noiseA
        xB_new = np.dot(self.A,self.xB[-1]) + noiseB
        
        self.xA = np.vstack((self.xA,[xA_new]))
        self.xB = np.vstack((self.xB,[xB_new]))
        
        self.__current_generation += 1
        
        return xA_new, xB_new

    def reset(self):
        if self.experimenttype == 'sisters':
            start   = self.random()
            self.xA = np.array([start])
            self.xB = np.array([start])
        else:
            self.xA = np.array([self.random()])
            self.xB = np.array([self.random()])
        
        self.__current_generation = 0


    def run(self,steps = None):
        self.reset()
        if steps is None:
            steps = self.default_steps
            
        for i in np.arange(steps):
            self.step()
        return self.xA, self.xB

        
    def rotation(self,angle):
        # angle in (0 ... 1)
        return np.array([[np.cos(2*np.pi*angle),np.sin(2*np.pi*angle)],[-np.sin(2*np.pi*angle),np.cos(2*np.pi*angle)]],dtype=np.float)


    def interval(self,value,minval = 0,maxval = 1):
        step = maxval - minval
        cval = value
        while cval < minval:cval += step
        while cval > maxval:cval -= step
        return cval


    def projection(self,x, alphaangle = None):
        if alphaangle is None:  alphavec = self.alpha
        else:                   alphavec = np.dot(self.rotation(alphaangle),self.__A_eigvec[0])
        
        alphavec = alphavec/np.linalg.norm(alphavec)
        return np.dot(alphavec,x)
    
    
    def output(self,step = None,alphaangle = None):
        if step is None:    step = -1
        else:               step = np.max([0,np.min([len(self.xA)-1,step])])
        return self.projection(self.xA[step],alphaangle), self.projection(self.xB[step],alphaangle)


    def output_all(self,alphaangle = None):
        return np.array([self.output(step = i,alphaangle = alphaangle) for i in range(len(self))]).T
    
    
    def __len__(self):
        assert len(self.xA) == len(self.xB)
        return len(self.xA)


    # analytical results
    def compute_AmATk(self,m = 0,k = 0):
        # exact result for product A^m (A.T)^k in our parameterization, using that first eigenvector is (0,1)
        l1m = np.power(self.__A_eigval[0],m)
        l2m = np.power(self.__A_eigval[1],m)
        l1k = np.power(self.__A_eigval[0],k)
        l2k = np.power(self.__A_eigval[1],k)
        
        return np.array([[ l2m * l2k,                                            l2m*(l1k-l2k)/np.tan(2 * np.pi * self.__A_angle) ],
                         [ l2k * (l1m - l2m)/np.tan(2 * np.pi * self.__A_angle), l1m*l1k + (l1k-l2k)*(l1m-l2m)/(np.tan(2 * np.pi * self.__A_angle)**2)]], dtype = np.float)

    def compute_Am(self,m=1):
        l1m = np.power(self.__A_eigval[0],m)
        l2m = np.power(self.__A_eigval[1],m)
        return np.array([[l2m,0],[(l1m-l2m)/np.tan(2 * np.pi * self.__A_angle),l1m]],dtype=np.float)
    
    def compute_ATk(self,k=1):
        l1k = np.power(self.__A_eigval[0],k)
        l2k = np.power(self.__A_eigval[1],k)
        return np.array([[l2k,(l1k-l2k)/np.tan(2 * np.pi * self.__A_angle)],[0,l1k]],dtype=np.float)
        

    def compute_sumInf_AmATm(self):
        # sum_{m=0}^Infinity A^m (A.T)^m
        il1l1 = 1./(1-self.__A_eigval[0]**2)
        il2l2 = 1./(1-self.__A_eigval[1]**2)
        il1l2 = 1./(1-self.__A_eigval[0]*self.__A_eigval[1])
        itan  = 1./np.tan(2 * np.pi * self.__A_angle)
        return np.array([[ il2l2,                  (il1l2 - il2l2) * itan],
                         [ (il1l2 - il2l2) * itan, il1l1 + (il1l1 - 2*il1l2 + il2l2)*itan]], dtype = np.float)
    

    def StationaryDifferenceCorrelations(self):
        if self.experimenttype == 'sisters':
            return np.zeros((2,2))
        elif self.experimenttype == 'nonsisters':
            return self.noiseamplitudes[0]**2 * self.sumInf_AmATm
        elif self.experimenttype == 'control':
            return (self.noiseamplitudes[0]**2 + self.noiseamplitudes[1]**2) * self.sumInf_AmATm


    def VarianceDifferenceProjectionValue(self,generation = 1):
        # Var[alpha.T sum(xA - xB) sum(xA - xB).T alpha]
        Adx0AT = np.matmul(self.A,np.matmul(self.StationaryDifferenceCorrelations(),self.A.T))
        noiseamplitude2 = self.noiseamplitudes[0]**2
        if self.experimenttype == 'control':
            noiseamplitude2 += self.noiseamplitudes[1]**2
        corr = np.zeros((2,2))
        for m in range(generation):
            for k in range(generation):
                corr += np.matmul(self.compute_Am(m),np.matmul(Adx0AT + noiseamplitude2 * (generation - np.max((m,k))) * np.eye(2),self.compute_ATk(k)))
        return np.dot(self.alpha,np.dot(corr,self.alpha))


    def VarianceDifferenceProjection(self,generation = 1):
        if len(self.vardiff) <= generation:
            tmp = list()
            for i in range(len(self.vardiff),generation + 1):
                tmp.append(self.VarianceDifferenceProjectionValue(generation = i))
            self.vardiff = np.concatenate([self.vardiff,tmp])
        return self.vardiff[:generation+1] # array starts at generation 0
        

        

def main():
    parser = argparse.ArgumentParser()
    parser_process = parser.add_argument_group(description = "==== process parameters ====")
    parser_process.add_argument("-n","--steps",default=10,type=int)
    parser_process.add_argument("-T","--experimenttype",choices = ('sisters','nonsisters','control'),default='sisters')
    
    parser_single = parser.add_argument_group(description = "==== single step parameters ====")
    parser_single.add_argument("-A","--A_angle",default=.33,type=float)
    parser_single.add_argument("-E","--eigval",type=float,nargs=2,default=[0.9,0.5])
    parser_single.add_argument("-a","--alpha_angle",default=.67,type=float)
    parser_single.add_argument("-N","--noiseamplitudes",default=[1,1],nargs=2,type=float)

    args = parser.parse_args()

    p = ARP(**vars(args))

    for s in np.arange(args.steps):
        print('{:4d} {:8.4f} {:8.4f}'.format(s,*p.output()))
        p.step()
        

if __name__ == "__main__":
    main()
