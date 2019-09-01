#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math
import pandas as pd


class MDARP(object):
    """
    multi-dimensional autoregressive process with several noise terms
    """
    def __init__(self,**kwargs):
        # construct matrix A from eigenvalues and angle (0...1) between eigenvectors, assume first direction is (0,1)
        self.__A_eigval    = np.array([self.interval(x) for x in kwargs.get('lambda',[0.9,0.5])],dtype=np.float)
        self.__A_angle     = self.interval(kwargs.get('A_angle',0.33)) * 2 * np.pi
        evec1              = np.array([0,1],dtype=np.float)
        evec2              = np.dot(self.rotation(self.__A_angle),evec1)
        self.__A_eigvec    = (evec1 / np.linalg.norm(evec1), evec2 / np.linalg.norm(evec2))
        evmat              = np.array(self.__A_eigvec,dtype=float).T
        self.A             = np.matmul(np.matmul(evmat,np.diag(self.__A_eigval)), np.linalg.inv(evmat))
        
        # projection vector is also with respect to first EVec
        self.__alpha_angle = self.interval(kwargs.get('alpha_angle',0.67)) * 2 * np.pi
        self.alpha         = np.dot(self.rotation(self.__alpha_angle),evec1/np.linalg.norm(evec1))

        # noise, linage, env
        self.noiseamplitudes = np.array(kwargs.get('noiseamplitudes',[1,1,1]),dtype=np.float)

        self.xA = self.random()
        self.xB = self.random()
        
        self.cumu_xA = np.zeros(2)
        self.cumu_xB = np.zeros(2)

        self.experimenttype = kwargs.get('experimenttype','sisters')
        
        self.__curgen = 0


    def random(self,mean = 0, var = 1):
        return np.random.normal(size=2)


    def step(self):
        noiseA = self.noiseamplitudes[0] * self.random()
        noiseB = self.noiseamplitudes[0] * self.random()
        if self.experimenttype == 'sisters' or self.experimenttype == 'nonsisters':
            if self.__curgen == 0 and self.experimenttype == 'sisters':
                xiL = self.random()
                noiseA += self.noiseamplitudes[1] * xiL
                noiseB += self.noiseamplitudes[1] * xiL
            else:
                noiseA += self.noiseamplitudes[1] * self.random()
                noiseB += self.noiseamplitudes[1] * self.random()
            xiE = self.random()
            
            noiseA += self.noiseamplitudes[2] * xiE
            noiseB += self.noiseamplitudes[2] * xiE
        else:
            noiseA += self.noiseamplitudes[1] * self.random() + self.noiseamplitudes[2] * self.random()
            noiseB += self.noiseamplitudes[1] * self.random() + self.noiseamplitudes[2] * self.random()

        self.xA = np.dot(self.A,self.xA) + noiseA
        self.xB = np.dot(self.A,self.xB) + noiseB
        
        self.cumu_xA += self.xA
        self.cumu_xB += self.xB
        
        self.__curgen += 1
        
        return self.xA,self.xB

    def run(self,steps):
        for i in np.arange(steps):
            self.step()
        return self.xA, self.xB

        
    def rotation(self,angle):
        return np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]],dtype=np.float)


    def interval(self,value,minval = 0,maxval = 1):
        step = maxval - minval
        cval = value
        while cval < minval:cval += step
        while cval > maxval:cval -= step
        return cval


    def projection(self,x):
        return np.dot(self.alpha,x)
    
    
    def output(self,cumulative = False):
        if cumulative:  return self.projection(self.cumu_xA), self.projection(self.cumu_xB)
        else:           return self.projection(self.xA), self.projection(self.xB)



def main():
    parser = argparse.ArgumentParser()
    parser_process = parser.add_argument_group(description = "==== process parameters ====")
    parser_process.add_argument("-n","--steps",default=10,type=int)
    parser_process.add_argument("-T","--experimenttype",choices = ('sisters','nonsisters','control'),default='sisters')
    
    parser_single = parser.add_argument_group(description = "==== single step parameters ====")
    parser_single.add_argument("-A","--A_angle",default=.33,type=float)
    parser_single.add_argument("-L","--lambda",type=float,nargs=2,default=[0.9,0.5])
    parser_single.add_argument("-a","--alpha_angle",default=.67,type=float)
    parser_single.add_argument("-N","--noiseamplitudes",default=[1,1,1],nargs=3,type=float)

    args = parser.parse_args()

    p = MDARP(**vars(args))

    for s in np.arange(args.steps):
        print('{:4d} {:8.4f} {:8.4f}'.format(s,*p.output()))
        p.step()
        

if __name__ == "__main__":
    main()
