#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import numpy as np

RATIO_MEAN=5
RATIO_SD=1

#def sampleBeta(k,m,mode,conc):
#    alpha=mode*(conc-2)+1
#    beta=(1-mode)*(conc-2)+1
#    s=np.random.beta(k+alpha,m+beta)
#    return s

class pool:
    def __init__(self,freq,dnaRef,dnaAlt,rnaRef,rnaAlt):
        self.freq=freq
        self.dnaRef=dnaRef
        self.dnaAlt=dnaAlt
        self.rnaRef=rnaRef
        self.rnaAlt=rnaAlt

def simGenotypes(n,popFreq):
    genotypes=[]
    for i in range(n):
        allele1=1 if np.random.uniform()<popFreq else 0
        allele2=1 if np.random.uniform()<popFreq else 0
        gt=[allele1,allele2]
        genotypes.append(gt)
    return genotypes

def hasBothAlleles(genotypes):
    (numRef,numAlt)=countAlleles(genotypes)
    return numRef>0 and numAlt>0

def countAlleles(genotypes):
    numRef=0; numAlt=0
    for gt in genotypes:
        if(gt[0]==0): numRef+=1
        else: numAlt+=1
        if(gt[1]==0): numRef+=1
        else: numAlt+=1
    return (numRef,numAlt)

def splitIntoPools(genotypes,NUM_POOLS,indivPerPool):
    pools=[]
    nextIndiv=0
    for i in range(NUM_POOLS):
        pool=genotypes[nextIndiv:(nextIndiv+indivPerPool)]
        nextIndiv+=indivPerPool
        pools.append(pool)
    return pools

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <#variants> <#indiv> <#pools> <TOTAL_COUNT>\n")
(NUM_VARIANTS,NUM_INDIV,NUM_POOLS,TOTAL_COUNT)=sys.argv[1:]
NUM_VARIANTS=int(NUM_VARIANTS)
NUM_INDIV=int(NUM_INDIV)
NUM_POOLS=int(NUM_POOLS)
if(NUM_INDIV % NUM_POOLS != 0):
    raise Exception("#indiv must be divisib by #pools")
indivPerPool=int(NUM_INDIV/NUM_POOLS)
TOTAL_COUNT=int(TOTAL_COUNT)

for i in range(NUM_VARIANTS):
    theta=np.random.lognormal(0,1)
    r_ref=np.random.normal(RATIO_MEAN,RATIO_SD)
    while(r_ref==0): r_ref=np.random.normal(RATIO_MEAN,RATIO_SD)
    r_alt=r_ref*theta
    popFreq=None; genotypes=None
    while(True):
        popFreq=np.random.uniform(0,0.5)
        genotypes=simGenotypes(NUM_INDIV,popFreq)
        if(hasBothAlleles(genotypes)): break
    pools=splitIntoPools(genotypes,NUM_POOLS,indivPerPool)
    for pool in pools:
        (numRef,numAlt)=countAlleles(pool)
        localFreq=float(numAlt)/float(numAlt+numRef)
        print(localFreq)
    
#(variant (id chr1:100038008:T:C)
#        (pool 1
#                (freq 0.4)
#                (DNA (ref 55) (alt 50))
#                (RNA (ref 261) (alt 227))
#        )



