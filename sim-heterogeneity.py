#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import random
import gzip
import statistics
import numpy as np
import ProgramName
from Rex import Rex

SEED=0
random.seed(SEED)
rng=np.random.default_rng(seed=SEED)
rex=Rex()

class Donor:
    def __init__(self):
        self.genotypes=[] # indexed by variant, each genotype=[int,int]
    def numVariants(self):
        return len(self.genotypes)

class Pool:
    def __init__(self):
        self.donors=[]
        self.AFs=None # indexed by variant
    def numVariants(self):
        if(len(self.donors)==0): return 0
        return self.donors[0].numVariants()
    def numDonors(self):
        return len(self.donors)
    def computeAFs(self):
        numVariants=self.numVariants()
        self.AFs=[0]*numVariants
        for i in range(numVariants):
            totalAlleles=0
            totalAlts=0
            for donor in self.donors:
                gt=donor.genotypes[i]
                totalAlleles+=2
                totalAlts+=gt[0]+gt[1]
            self.AFs[i]=float(totalAlts)/float(totalAlleles)

def computeHeterogeneity(pools):
    values=[pool.AFs[0] for pool in pools]
    return statistics.variance(values) if len(values)>1 else 0
    
def computePoolAFs(pools):
    if(len(pools)==0): raise Exception("no pools")
    numVariants=pools[0].numVariants()
    for pool in pools: pool.computeAFs()
        
def initPools(numPools):
    return [Pool() for _ in range(numPools)]
        
def initDonors(numDonors):
    return [Donor() for _ in range(numDonors)]
        
def loadVCF(vcfFilename,maxVariants,minAF,maxAF):
    numLoaded=0
    donors=None
    with gzip.open(vcfFilename,"rt") as IN:
        # For each variant:
        for line in IN:
            if(line[0]=="#"): continue
            fields=line.rstrip().split()
            numFields=len(fields)
            if(numFields<100): continue
            rex.findOrDie("AF=([^;]+);",fields[7])
            AF=float(rex[1])
            if(AF<minAF or AF>maxAF): continue # allele frequency filter
            numDonors=numFields-9
            if(donors is None): donors=initDonors(numDonors)
            # Append genotypes of this variant to each donor
            for i in range(9,numFields):
                gt=fields[i]
                rex.findOrDie("(\d)\|(\d)",gt)
                gt=[int(rex[1]),int(rex[2])]
                donors[i-9].genotypes.append(gt)
            numLoaded+=1
            if(numLoaded>=maxVariants): break
    return donors

def assignDonorsToPools(donors,pools):
    numPools=len(pools)
    for index, element in enumerate(donors):
        pools[index % numPools].donors.append(element)

def addNoise(x,sd):
    if(x<sd): sd=x/2
    noise=rng.normal(loc=0,scale=sd)
    newValue=x+noise
    while(newValue<=0):
        sd/=2
        if(sd<0.01): return x
        noise=rng.normal(loc=0,scale=sd)
        newValue=x+noise
    if(newValue>1): newValue=1
    newValue=round(newValue,4)
    return newValue
        
def simCounts(numVariants,pools,readsPerVariant,theta,numReps,poolsOrReps,
              sd,heterogeneity):
    multiplePools=(poolsOrReps=="pools")
    numPools=len(pools)
    numGroups=numPools if poolsOrReps=="pools" else numReps
    print("variants=",numVariants," groups=",numGroups,
          " heterogeneity=",heterogeneity," theta=",theta,
          " each group is (ref,alt,AF)",sep="")
    for i in range(numVariants):
        print("var",i+1,"\t",sep="",end="")
        for j in range(numPools):
            pool=pools[j]
            for k in range(numReps):
                varIndex=0 # Always simulate from the same real variant!
                generateGroup(pool,readsPerVariant,theta,varIndex,j,sd)
                if(multiplePools and j<numPools-1 or
                   not multiplePools and k<numReps-1):
                    print("\t",end="")
        print()

def generateGroup(pool,readsPerVariant,theta,varIndex,poolIndex,sd):
    p=pool.AFs[varIndex]
    p=addNoise(p,sd) # only added to DNA: represents transfection drift
    alt=rng.binomial(readsPerVariant,p,1)[0]
    ref=readsPerVariant-alt
    print("dna=",end="")
    print(ref,",",alt,",",p,"\t",sep="",end="")
    print("rna=",end="")
    q=theta*p/(1-p+theta*p)
    q=round(q,5)
    alt=rng.binomial(readsPerVariant,q,1)[0]
    ref=readsPerVariant-alt
    print(ref,",",alt,",",q,"\t",sep="",end="")
            
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=10):
    exit(ProgramName.get()+" <in.vcf.gz> <num-variants> <min-AF> <max-AF> <num-pools-or-reps> <pools|reps> <theta> <#reads-per-variant-per-pool> <noiseSD>\n")
(vcfFilename,numVariants,minAF,maxAF,numGroups,poolsOrReps,theta,readsPerVariant,noiseSD)=sys.argv[1:]

# Parse command line
numVariants=int(numVariants)
minAF=float(minAF)
maxAF=float(maxAF)
numGroups=int(numGroups)
theta=float(theta)
readsPerVariant=int(readsPerVariant)
noiseSD=float(noiseSD)
numPools=1
numReps=1
if(poolsOrReps=="pools"): numPools=numGroups
elif(poolsOrReps=="reps"): numReps=numGroups
else: raise Exception("specify pools or reps")

# Load VCF file
MAX_VARIANTS=1 # always simulate from the same real variant!
donors=loadVCF(vcfFilename,MAX_VARIANTS,minAF,maxAF)
numDonors=len(donors)

# Assign donors to pools
pools=initPools(numPools)
assignDonorsToPools(donors,pools)
computePoolAFs(pools)
heterogeneity=computeHeterogeneity(pools)

# Simulate counts
simCounts(numVariants,pools,readsPerVariant,theta,numReps,poolsOrReps,
          noiseSD,heterogeneity)
