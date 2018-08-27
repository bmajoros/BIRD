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

def writeReadCounts(fields,start,numReps,varName,OUT):
    print(varName,"<- c(",file=OUT,end="")
    for rep in range(numReps):
        print(fields[start+rep*2],file=OUT,end="")
        if(rep+1<numReps): print(",",file=OUT,end="")
    print(")",file=OUT)

def writeInitializationFile(fields,filename):
    DNAreps=int(fields[1])
    totalRef=0; totalAlt=0
    for i in range(DNAreps):
        totalRef+=int(fields[2+i*2])
        totalAlt+=int(fields[3+i*2])
    v=float(totalAlt)/float(totalAlt+totalRef)
    if(v==0): v=0.01
    rnaIndex=2+2*DNAreps
    RNAreps=int(fields[rnaIndex])
    OUT=open(filename,"wt")
    print("p <-",v,file=OUT)
    print("theta <- 1",file=OUT)
    print("qi <- c(",file=OUT,end="")
    for i in range(RNAreps-1):
        print(v,",",sep="",end="",file=OUT)
    print(v,")",sep="",file=OUT)
    OUT.close()

def writeInputsFile(fields,filename):
    DNAreps=int(fields[1])
    rnaIndex=2+2*DNAreps
    RNAreps=int(fields[rnaIndex])
    OUT=open(filename,"wt")
    print("N_DNA <-",str(DNAreps),file=OUT)
    writeReadCounts(fields,3,DNAreps,"a",OUT) # alt
    writeReadCounts(fields,2,DNAreps,"b",OUT) # ref
    print("N_RNA <-",str(RNAreps),file=OUT)
    writeReadCounts(fields,rnaIndex+2,RNAreps,"k",OUT) # alt
    writeReadCounts(fields,rnaIndex+1,RNAreps,"m",OUT) # ref
    OUT.close()

def getMedian(thetas):
    # Precondition: thetas is already sorted
    n=len(thetas)
    mid=int(n/2)
    if(n%2==0): return (thetas[mid-1]+thetas[mid])/2.0
    return thetas[mid]

def getCredibleInterval(thetas,alpha):
    halfAlpha=alpha/2.0
    n=len(thetas)
    leftIndex=int(halfAlpha*n)
    rightIndex=n-leftIndex
    left=thetas[leftIndex+1]
    right=thetas[rightIndex-1]
    return (left,right)

def runVariant(model,fields,numSamples,outfile):
    # Write inputs file for STAN
    if(len(fields)<10): return (None,None)
    writeInputsFile(fields,INPUT_FILE)
    writeInitializationFile(fields,INIT_FILE)

    # Run STAN model
    init=" init="+INIT_FILE
    cmd=model+" sample thin=1"+\
        " num_samples="+numSamples+\
        " num_warmup="+str(WARMUP)+\
        " data file="+INPUT_FILE+\
        init+\
        " output file="+OUTPUT_TEMP+" refresh=0 > "+STDERR
    os.system(cmd)

    # Parse MCMC output
    thetas=[];
    OUT=open(outfile,"wt")
    with open(OUTPUT_TEMP,"rt") as IN:
        for line in IN:
            if(len(line)==0 or line[0]=="#"): continue
            fields=line.rstrip().split(",")
            numFields=len(fields)
            if(numFields>0 and fields[0]=="lp__"):
                printFields(fields,OUT)
                thetaIndex=getFieldIndex("theta",fields)
            else:
                writeToFile(fields,OUT)
                theta=float(fields[thetaIndex])
                thetas.append(theta)
    OUT.close()
    thetas.sort(key=lambda x: x)
    return thetas

def summarize(thetas,fields,ID,minRight):
    maxLeft=1.0/minRight
    n=len(thetas)
    median=getMedian(thetas)
    (CI_left,CI_right)=getCredibleInterval(thetas,ALPHA)
    reduction=0
    increase=0
    for i in range(n):
        if(thetas[i]<maxLeft): reduction+=1
        if(thetas[i]>minRight): increase+=1
    leftP=float(reduction)/float(n)
    rightP=float(increase)/float(n)
    Preg=leftP if leftP>rightP else rightP
    print(ID,median,CI_left,CI_right,Preg,sep="\t")

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:")
if(len(args)!=6):
    exit(ProgramName.get()+" [-s file] <model> <min-effect> <input.txt> <output.txt> <#MCMC-samples> <firstVariant-lastVariant>\n   -s = save raw STAN file\n   variant range is zero-based and inclusive\n   min-effect (lambda) must be >= 1\n")
(model,minEffect,inFile,outfile,numSamples,numVariants)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
if(not rex.find("(\d+)-(\d+)",numVariants)):
    exit(numVariants+": specify range of variants: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])
minEffect=float(minEffect)
if(minEffect<1): raise Exception("Min-effect must be >= 1")

# Process all input lines, each line = one variant (one MCMC run)
thetaIndex=None
variantIndex=0
with open(inFile,"rt") as IN:
    for line in IN:
        # Check whether this variant is in the range to be processed
        if(variantIndex<firstIndex):
            variantIndex+=1
            continue
        elif(variantIndex>lastIndex): break
        fields=line.rstrip().split()
        ID=fields[0]
        thetas=runVariant(model,fields,numSamples,outfile)
        if(thetas is None): continue
        summarize(thetas,fields,ID,minEffect)
        variantIndex+=1
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)

