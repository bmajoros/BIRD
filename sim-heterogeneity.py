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

def initDonors(numDonors):
    donors=[]
    for i in range(numDonors):
        donors.append(Donor())
    return donors
        
def loadVCF(vcfFilename,maxVariants,minAF,maxAF):
    numLoaded=0
    donors=None
    with gzip.open(vcfFilename,"rt") as IN:
        for line in IN:
            if(line[0]=="#"): continue
            fields=line.rstrip().split()
            numFields=len(fields)
            if(numFields<100): continue
            rex.findOrDie("AF=([^;]+);",fields[7])
            AF=float(rex[1])
            if(AF<minAF or AF>maxAF): continue
            numLoaded+=1
            if(numLoaded>=maxVariants): break
            numDonors=numFields-9
            if(donors is None): donors=initDonors(numDonors)
            for i in range(9,numFields):
                gt=fields[i]
                rex.findOrDie("(\d)\|(\d)",gt)
                gt=[rex[1],rex[2]]
                donors[i-9].genotypes.append(gt)
    return donors
                    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <in.vcf.gz> <max-variants> <min-AF> <max-AF>\n")
(vcfFilename,maxVariants,minAF,maxAF)=sys.argv[1:]

maxVariants=int(maxVariants)
minAF=float(minAF)
maxAF=float(maxAF)
donors=loadVCF(vcfFilename,maxVariants,minAF,maxAF)
print(donors[0].genotypes)

