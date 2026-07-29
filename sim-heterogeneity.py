#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import random
import gzip
import ProgramName
from Rex import Rex

random.seed(0)
rex=Rex()

class Donor:
    def __init__(self):
        self.genotypes=[]

class Pool:
    def __init__(self):
        self.donors=[]

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
                gt=[rex[1],rex[2]]
                donors[i-9].genotypes.append(gt)
    return donors

def assignDonorsToPools(donors,pools):
    numPools=len(pools)
    for index, element in enumerate(donors):
        pools[index % numPools].donors.append(element)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <in.vcf.gz> <max-variants> <min-AF> <max-AF> <num-pools>\n")
(vcfFilename,maxVariants,minAF,maxAF,numPools)=sys.argv[1:]

# Parse command line
maxVariants=int(maxVariants)
minAF=float(minAF)
maxAF=float(maxAF)
numPools=int(numPools)

# Load VCF file
donors=loadVCF(vcfFilename,maxVariants,minAF,maxAF)
numDonors=len(donors)
#print(donors[0].genotypes)

# Assign donors to pools
pools=initPools(numPools)
assignDonorsToPools(donors,pools)
#for pool in pools: print(len(pool.donors))

# Simulate counts
