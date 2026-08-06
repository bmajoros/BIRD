#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import gzip
import math
import statistics
from scipy.stats import binom
import ProgramName
from Rex import Rex

DEBUG=0
rex=Rex()

def processVCF(filename,MAX_VARIANTS,NUM_DONORS,MAX_POOLS,COVERAGE,
               MIN_AF,MAX_AF):
    numProcessed=0
    sims=[]
    with gzip.open(filename,"rt") as IN:
        for line in IN:
            if(len(line)<1): continue
            if(line[0]=='#'): continue
            genotypes=processLine(line,NUM_DONORS)
            if(genotypes is None): continue
            AF=getAF(genotypes)
            if(AF==0.0): continue
            if(AF<MIN_AF or AF>MAX_AF): continue
            #print("AF=",AF)
            dropoutProbs=getDropoutProbs(genotypes,NUM_DONORS,MAX_POOLS,
                                         COVERAGE)
            sims.append(dropoutProbs)
            numProcessed+=1
            if(numProcessed>=MAX_VARIANTS): break
    print("NumPools\tDropout")
    for numPools in range(MAX_POOLS):
        ave=statistics.mean([probs[numPools] for probs in sims])
        print(numPools+1,round(ave,5),sep="\t")

def getAF(genotypes):
    A=countAlts(genotypes)
    N=2*len(genotypes)
    return float(A)/float(N)
        
def getDropoutProbs(genotypes,NUM_DONORS,MAX_POOLS,COVERAGE):
    probs=[]
    for numPools in range(1,MAX_POOLS+1):
        partitioning=partition(genotypes,numPools)
        p=getDropout(partitioning,COVERAGE)
        probs.append(p)
    return probs

def countAlts(pool):
    return sum([gt[0]+gt[1] for gt in pool])

def getDropout(pools,COVERAGE):
    if(DEBUG): print("NPOOLS=",len(pools))
    product=1
    for pool in pools:
        N=2*len(pool)
        freq=getAF(pool)
        p=binom.pmf(0,COVERAGE,freq)
        #print("binom(0,",N,",",freq,")=",p)
        if(DEBUG):
            print("\tAF=",freq,"A=",countAlts(pool),"N=",N,sep="\t",end="")
            print("\tbinom=",p)
        if(math.isnan(p)): raise Exception("NaN in binom(",0,N,freq)
        product*=p
    if(DEBUG): print("\tprod=",product)
    return product

def partition(genotypes,NUM_POOLS):
    pools=[[] for _ in range(NUM_POOLS)]
    for i, gt in enumerate(genotypes):
        pools[i % NUM_POOLS].append(gt)
    return pools
            
def processLine(line,NUM_DONORS):
    fields=line.rstrip().split()
    if(len(fields)<9+NUM_DONORS): return None
    genotypes=fields[9:(9+NUM_DONORS)]
    GT=[]
    for genotype in genotypes:
        rex.findOrDie("(\d)\|(\d)",genotype)
        gt=[int(rex[1]),int(rex[2])]
        if(gt[0] not in [0,1] or gt[1] not in [0,1]): return None
        GT.append(gt)
    return GT
            
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <in.vcf.gz> <max-pools> <max-variants> <#donors> <coverage> <min-AF> <max-AF>\nwhere coverage is the average # fragments per position\n")
(vcfFilename,MAX_POOLS,MAX_VARIANTS,NUM_DONORS,COVERAGE,
 MIN_AF,MAX_AF)=sys.argv[1:]
MAX_POOLS=int(MAX_POOLS)
MAX_VARIANTS=int(MAX_VARIANTS)
NUM_DONORS=int(NUM_DONORS)
COVERAGE=int(COVERAGE)
MIN_AF=float(MIN_AF)
MAX_AF=float(MAX_AF)

processVCF(vcfFilename,MAX_VARIANTS,NUM_DONORS,MAX_POOLS,COVERAGE,
           MIN_AF,MAX_AF)



