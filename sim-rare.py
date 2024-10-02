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

MAX_FREQ=0.01
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

def sampleBinomial(n,p):
    altCount=np.random.binomial(n,p)
    refCount=n-altCount
    return (refCount,altCount)

def sampleNB(TOTAL_DNA,ratio):
    ALPHA=1 # 1e-10
    BETA=0.05 # 1e-10
    p=(BETA+1)/(BETA+ratio+1)
    r=TOTAL_DNA+ALPHA
    TOTAL_RNA=np.random.negative_binomial(r,p)
    return TOTAL_RNA

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=7):
    exit(ProgramName.get()+" <#variants> <#indiv> <#pools> <TOTAL_COUNT> <out:data> <out:truth>\n")
(NUM_VARIANTS,NUM_INDIV,NUM_POOLS,TOTAL_COUNT,outDataFile,outTruthFile)=\
    sys.argv[1:]
NUM_VARIANTS=int(NUM_VARIANTS)
NUM_INDIV=int(NUM_INDIV)
NUM_POOLS=int(NUM_POOLS)
if(NUM_INDIV % NUM_POOLS != 0):
    raise Exception("#indiv must be divisib by #pools")
indivPerPool=int(NUM_INDIV/NUM_POOLS)
TOTAL_DNA=int(TOTAL_COUNT)

OUT_DATA=open(outDataFile,"wt")
OUT_TRUTH=open(outTruthFile,"wt")
for i in range(NUM_VARIANTS):
    variantID="variant"+str(i+1)
    print("(variant (id ",variantID,")",sep="",file=OUT_DATA)
    theta=np.random.lognormal(0,1)
    r_ref=np.random.normal(RATIO_MEAN,RATIO_SD)
    while(r_ref==0): r_ref=np.random.normal(RATIO_MEAN,RATIO_SD)
    r_alt=r_ref*theta
    print(variantID,round(theta,3),
          round(r_ref,3),round(r_alt,3),
          sep="\t",file=OUT_TRUTH)
    #print("\t(theta ",round(theta,3),")",sep="",file=OUT_DATA)
    #print("\t(r_ref ",round(r_ref,3),")",sep="",file=OUT_DATA)
    #print("\t(r_alt ",round(r_alt,3),")",sep="",file=OUT_DATA)
    popFreq=None; genotypes=None
    while(True):
        popFreq=np.random.uniform(0,MAX_FREQ)
        genotypes=simGenotypes(NUM_INDIV,popFreq)
        if(hasBothAlleles(genotypes)): break
    pools=splitIntoPools(genotypes,NUM_POOLS,indivPerPool)
    poolID=0
    for pool in pools:
        poolID+=1
        (numRef,numAlt)=countAlleles(pool)
        localFreq=float(numAlt)/float(numAlt+numRef)
        print("\t(pool ",poolID," (freq ",localFreq,")",sep="",file=OUT_DATA)
        TOTAL_RNA=int(TOTAL_DNA*(r_ref*(1-localFreq) +r_alt*localFreq))
        dnaRef=None;dnaAlt=None;rnaRef=None;rnaAlt=None
        if(localFreq>0 and localFreq<1): # Het pool
            p=localFreq
            q=theta*p/(1-p+theta*p)
            (dnaRef,dnaAlt)=sampleBinomial(TOTAL_DNA,p)
            (rnaRef,rnaAlt)=sampleBinomial(TOTAL_RNA,q)
        elif(localFreq==0): # Homozygous reference
            TOTAL_RNA=sampleNB(TOTAL_DNA,r_ref)
            dnaRef=TOTAL_DNA; dnaAlt=0
            rnaRef=TOTAL_RNA; rnaAlt=0
        else: # Homozygous alternate
            TOTAL_RNA=sampleNB(TOTAL_DNA,r_alt)
            dnaAlt=TOTAL_DNA; dnaRef=0
            rnaAlt=TOTAL_RNA; rnaRef=0
        print("\t\t(DNA (ref ",dnaRef,") (alt ",dnaAlt,"))",sep="",
              file=OUT_DATA)
        print("\t\t(RNA (ref ",rnaRef,") (alt ",rnaAlt,"))",sep="",
              file=OUT_DATA)
        print("\t)",file=OUT_DATA)
    print(")",file=OUT_DATA)
    
#(variant (id chr1:100038008:T:C)
#        (pool 1
#                (freq 0.4)
#                (DNA (ref 55) (alt 50))
#                (RNA (ref 261) (alt 227))
#        )



