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

def fixP(p):
    newP=[x if x>0 and x<1 else 0.5 for x in p]
    return newP
    
def writeInitializationFile(stan,variant,filename):
    OUT=open(filename,"wt")
    print("theta <- 1",file=OUT)
    print("r_ref <- 1",file=OUT)
    print("s <- 1",file=OUT)
    freqs=variant.getFreqs()
    freqs=fixP(freqs)
    numPools=variant.numPools()
    stan.writeOneDimArray("p",freqs,numPools,OUT)
    OUT.close()

def writeInputsFile(stan,variant,filename):
    OUT=open(filename,"wt")
    numPools=variant.numPools()
    print("N_POOLS <-",numPools,file=OUT)
    poolTypes=[pool.getPoolType() for pool in variant.pools]
    stan.writeOneDimArray("POOL_TYPE",poolTypes,len(poolTypes),OUT)
    freqs=variant.getFreqs()
    stan.writeOneDimArray("pop_freq",freqs,numPools,OUT)
    print("pop_conc <- ",POP_CONC,file=OUT)
    dnaAltCounts=[pool.DNA[0].alt for pool in variant.pools]
    dnaRefCounts=[pool.DNA[0].ref for pool in variant.pools]
    rnaAltCounts=[pool.RNA[0].alt for pool in variant.pools]
    rnaRefCounts=[pool.RNA[0].ref for pool in variant.pools]
    stan.writeOneDimArray("a",dnaAltCounts,numPools,OUT)
    stan.writeOneDimArray("b",dnaRefCounts,numPools,OUT)
    stan.writeOneDimArray("k",rnaAltCounts,numPools,OUT)
    stan.writeOneDimArray("m",rnaRefCounts,numPools,OUT)
    print("mu <-",MU,file=OUT)
    print("sigma2 <-",SIGMA2,file=OUT)
    print("alpha <-",ALPHA,file=OUT)
    print("beta <-",BETA,file=OUT)
    OUT.close()

def runVariant(stan,variant,numSamples,outfile):
    # Write inputs file for STAN
    writeInputsFile(stan,variant,INPUT_FILE)
    writeInitializationFile(stan,variant,INIT_FILE)

    # Run STAN model
    cmd=stan.getCmd(WARMUP,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE)
    if(DEBUG):
        print(cmd)
        exit()
    os.system(cmd)

    # Parse MCMC output
    parser=StanParser(OUTPUT_TEMP)
    thetas=parser.getVariable("theta")    
    return (thetas,parser)

def summarize(parser,thetas,ID,minRight):
    (median,CI_left,CI_right)=parser.getMedianAndCI(0.95,"theta")
    maxLeft=1.0/minRight
    leftP=parser.getLeftTail("theta",maxLeft)
    rightP=parser.getRightTail("theta",minRight)
    Preg=leftP if leftP>rightP else rightP
    medianRatio=parser.getMedianAndCI(0.95,"r_ref")[0]
    #print(ID,round(median,3),round(CI_left,3),round(CI_right,3),
    #      round(Preg,3),sep="\t")
    print(ID,round(median,3),round(CI_left,3),round(CI_right,3),
          round(Preg,3),round(medianRatio,3),sep="\t")

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:t:")
if(len(args)!=7):
    exit(ProgramName.get()+" [-s stanfile] [-t thetafile] <model> <min-effect> <beta-concentration-parm> <input.essex> <output.txt> <#MCMC-samples> <firstVariant-lastVariant>\n   -s = save raw STAN file\n   -t = save theta samples\n   variant range is zero-based and inclusive\n   min-effect (lambda) must be >= 1\n")
(model,minEffect,POP_CONC,inFile,outfile,numSamples,numVariants)=args
stanFile=None
thetaFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
    if(key=="-t"): thetaFile=value
if(not rex.find(r"(\d+)-(\d+)",numVariants)):
    exit(numVariants+": specify range of variants: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])
minEffect=float(minEffect)
POP_CONC=float(POP_CONC)
if(minEffect<1): raise Exception("Min-effect must be >= 1")
THETA=None
if(thetaFile is not None): THETA=open(thetaFile,"wt")
stan=Stan(model)

# Process all input lines, each line = one variant (one MCMC run)
thetaIndex=None
variantIndex=0
pooledParser=PooledParser(inFile)
while(True):
    variant=pooledParser.nextVariant()
    if(variant is None): break
    # Check whether this variant is in the range to be processed
    if(variantIndex<firstIndex):
        variantIndex+=1
        continue
    elif(variantIndex>lastIndex): break
    #keep=variant.dropHomozygousPools()
    #if(not keep): continue
    (thetas,stanParser)=runVariant(stan,variant,numSamples,outfile)
    if(thetas is None): continue
    summarize(stanParser,thetas,variant.ID,minEffect)
    variantIndex+=1
    if(THETA is not None):
        for i in range(len(thetas)):
            print(thetas[i],file=THETA,end="")
            if(i<len(thetas)): print("\t",file=THETA,end="")
        print(file=THETA)
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)
if(THETA is not None): THETA.close()

