#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import random
import gzip
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
            numLoaded+=1
            if(numLoaded>=maxVariants): break
            numDonors=numFields-9
            if(donors is None): donors=initDonors(numDonors)
            # Append genotypes of this variant to each donor
            for i in range(9,numFields):
                gt=fields[i]
                rex.findOrDie("(\d)\|(\d)",gt)
                gt=[int(rex[1]),int(rex[2])]
                donors[i-9].genotypes.append(gt)
    return donors

def assignDonorsToPools(donors,pools):
    numPools=len(pools)
    for index, element in enumerate(donors):
        pools[index % numPools].donors.append(element)

def simCounts(numVariants,pools,readsPerVariant):
    numPools=len(pools)
    for i in range(numVariants):
        for j in range(numPools):
            pool=pools[j]
            AF=pool.AFs[i]
            alt=rng.binomial(readsPerVariant,AF,1)[0]
            ref=readsPerVariant-alt
            print(alt,",",ref)
            
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <in.vcf.gz> <num-variants> <min-AF> <max-AF> <num-pools> <theta> <#reads-per-variant-per-pool>\n")
(vcfFilename,maxVariants,minAF,maxAF,numPools,theta,readsPerVariant)=sys.argv[1:]

# Parse command line
maxVariants=int(maxVariants)
minAF=float(minAF)
maxAF=float(maxAF)
numPools=int(numPools)
theta=float(theta)
readsPerVariant=int(readsPerVariant)

# Load VCF file
donors=loadVCF(vcfFilename,maxVariants,minAF,maxAF)
numDonors=len(donors)
#print(donors[0].genotypes)

# Assign donors to pools
pools=initPools(numPools)
assignDonorsToPools(donors,pools)
computePoolAFs(pools)
#for pool in pools: print(len(pool.donors))
#print(pools[0].AFs)

# Simulate counts
simCounts(maxVariants,pools,readsPerVariant)
