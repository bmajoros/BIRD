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
ALPHA=0.05
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

def writeInitializationFile(stan,variant,filename):
    OUT=open(filename,"wt")
    #print("theta <- 1",file=OUT)
    freqs=variant.getFreqs()
    numPools=variant.numPools()
    stan.writeOneDimArray("p",freqs,numPools,OUT)
    #maxRnaReps=variant.getMaxRnaReps()
    #qiInit=[freqs]*maxRnaReps
    #stan.writeTwoDimArray("qi",qiInit,numPools,maxRnaReps,OUT)
    #print("c <- 100",file=OUT)
    #print("s <- 1",file=OUT)
    OUT.close()

def writeInputsFile(stan,variant,filename):
    OUT=open(filename,"wt")
    numPools=variant.numPools()
    print("N_POOLS <-",numPools,file=OUT)
    maxDnaReps=variant.getMaxDnaReps()
    #maxRnaReps=variant.getMaxRnaReps()
    print("MAX_DNA <- ",maxDnaReps,file=OUT)
    #print("MAX_RNA <- ",maxRnaReps,file=OUT)
    freqs=variant.getFreqs()
    stan.writeOneDimArray("pop_freq",freqs,numPools,OUT)
    dnaReps=variant.getDnaReps()
    #rnaReps=variant.getRnaReps()
    stan.writeOneDimArray("N_DNA",dnaReps,numPools,OUT)
    #stan.writeOneDimArray("N_RNA",rnaReps,numPools,OUT)
    dnaAltCounts=[[rep.alt for rep in pool.DNA] for pool in variant.pools]
    dnaRefCounts=[[rep.ref for rep in pool.DNA] for pool in variant.pools]
    #rnaAltCounts=[[rep.alt for rep in pool.RNA] for pool in variant.pools]
    #rnaRefCounts=[[rep.ref for rep in pool.RNA] for pool in variant.pools]
    dnaAltCounts=expandArray(dnaAltCounts,maxDnaReps)
    dnaRefCounts=expandArray(dnaRefCounts,maxDnaReps)
    #rnaAltCounts=expandArray(rnaAltCounts,maxRnaReps)
    #rnaRefCounts=expandArray(rnaRefCounts,maxRnaReps)
    stan.writeTwoDimArray("a",dnaAltCounts,numPools,maxDnaReps,OUT)
    stan.writeTwoDimArray("b",dnaRefCounts,numPools,maxDnaReps,OUT)
    #stan.writeTwoDimArray("k",rnaAltCounts,numPools,maxRnaReps,OUT)
    #stan.writeTwoDimArray("m",rnaRefCounts,numPools,maxRnaReps,OUT)
    OUT.close()

def expandArray(array,desiredSize):
    newArray=[]
    numPools=len(array)
    for i in range(numPools):
        oldPool=array[i]
        newPool=[x for x in oldPool]
        oldPoolSize=len(oldPool)
        for j in range(oldPoolSize,desiredSize):
            newPool.append(0)
        newArray.append(newPool)
    return newArray
    
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
    conc=parser.getVariable("pop_conc")    
    return (conc,parser)

#def summarize(parser,conc,ID,minRight):
#    (median,CI_left,CI_right)=parser.getMedianAndCI(1.0-ALPHA,"conc")
#    maxLeft=1.0/minRight
#    leftP=parser.getLeftTail("theta",maxLeft)
#    rightP=parser.getRightTail("theta",minRight)
#    Preg=leftP if leftP>rightP else rightP
#    print(ID,median,CI_left,CI_right,Preg,sep="\t")

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:t:")
if(len(args)!=6):
    exit(ProgramName.get()+" [-s stanfile] [-t thetafile] <model> <min-effect> <input.essex> <output.txt> <#MCMC-samples> <firstVariant-lastVariant>\n   -s = save raw STAN file\n   -t = save theta samples\n   variant range is zero-based and inclusive\n   min-effect (lambda) must be >= 1\n")
(model,minEffect,inFile,outfile,numSamples,numVariants)=args
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
    (conc,stanParser)=runVariant(stan,variant,numSamples,outfile)
    if(conc is None): continue
    #summarize(stanParser,thetas,variant.ID,minEffect)
    (median,CI_left,CI_right)=stanParser.getMedianAndCI(1.0-ALPHA,"pop_conc")
    print(variant.ID,median,sep="\t")
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

