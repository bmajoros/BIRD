#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2018 William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import math
import ProgramName
from Rex import Rex
rex=Rex()
import TempFilename
import getopt
from StanParser import StanParser
from PooledParser import PooledParser
from Stan import Stan

DEBUG=False
WARMUP=1000
STDERR=TempFilename.generate(".stderr")
INPUT_FILE=TempFilename.generate(".staninputs")
INIT_FILE=TempFilename.generate(".staninit")
OUTPUT_TEMP=TempFilename.generate(".stanoutputs")

def printFields(fields,hFile):
    numFields=len(fields)
    for i in range(7,numFields):
        print(i-6,"=",fields[i],sep="",end="",file=hFile)
        if(i<numFields-1): print("\t",end="",file=hFile)
    print(file=hFile)

def getFieldIndex(label,fields):
    numFields=len(fields)
    index=None
    for i in range(7,numFields):
        if(fields[i]==label): index=i
    return index

def writeToFile(fields,OUT):
    numFields=len(fields)
    for i in range(7,numFields):
        print(fields[i],end="",file=OUT)
        if(i<numFields-1): print("\t",end="",file=OUT)
    print(file=OUT)

def writeInitializationFile(stan,numVariants,filename):
    OUT=open(filename,"wt")
    print("mu <- 1",file=OUT)
    print("sigma2 <- 1",file=OUT)
    print("alpha <- 1",file=OUT)
    print("beta <- 1",file=OUT)
    initArray=[1]*numVariants
    stan.writeOneDimArray("r_ref",initArray,numVariants,OUT)
    stan.writeOneDimArray("lambda",initArray,numVariants,OUT)
    OUT.close()
    
def writeInputsFile(stan,b,m,filename):
    OUT=open(filename,"wt")
    numVariants=len(b)
    print("N_VARIANTS <-",numVariants,file=OUT)
    stan.writeOneDimArray("b",b,numVariants,OUT)
    stan.writeOneDimArray("m",m,numVariants,OUT)
    OUT.close()

def run(stan,b,m,numSamples):
    # Write inputs file for STAN
    writeInputsFile(stan,b,m,INPUT_FILE)
    writeInitializationFile(stan,len(b),INIT_FILE)

    # Run STAN model
    cmd=stan.getCmd(WARMUP,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE)
    if(DEBUG):
        print(cmd)
        exit()
    os.system(cmd)

    # Parse MCMC output
    parser=StanParser(OUTPUT_TEMP)
    return parser

def summarize(parser):
    (mu,mean,SD,min,max)=parser.getSummary("mu")
    (sigma2,mean,SD,min,max)=parser.getSummary("sigma2")
    (alpha,mean,SD,min,max)=parser.getSummary("alpha")
    (beta,mean,SD,min,max)=parser.getSummary("beta")
    print("mu=",round(mu,3),"sigma2=",round(sigma2,3),
          "alpha=",round(alpha,3),"beta=",round(beta,3),sep="\t")

def appendCounts(variant,dnaRefList,rnaRefList):
    newVar=variant.collapse()
    pool=newVar.pools[0]
    totalDnaRef=sum([rep.ref for rep in pool.DNA])
    totalRnaRef=sum([rep.ref for rep in pool.RNA])
    dnaRefList.append(totalDnaRef)
    rnaRefList.append(totalRnaRef)
    
#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:t:")
if(len(args)!=4):
    exit(ProgramName.get()+" [-s stanfile] [-t thetafile] <model> <input.essex> <#MCMC-samples> <firstVariant-lastVariant>\n   -s = save raw STAN file\n   variant range is zero-based and inclusive\n")
(model,inFile,numSamples,numVariants)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
if(not rex.find(r"(\d+)-(\d+)",numVariants)):
    exit(numVariants+": specify range of variants: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])
stan=Stan(model)

# Process all input lines, each line = one variant
variantIndex=0
pooledParser=PooledParser(inFile)
b=[]; m=[]
while(True):
    variant=pooledParser.nextVariant()
    if(variant is None): break
    # Check whether this variant is in the range to be processed
    if(variantIndex<firstIndex):
        variantIndex+=1
        continue
    elif(variantIndex>lastIndex): break
    appendCounts(variant,b,m)
    variantIndex+=1

# Run model on all variants jointly
stanParser=run(stan,b,m,numSamples)
summarize(stanParser)
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)

