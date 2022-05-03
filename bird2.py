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
from Stan import Stan
from StanParser import StanParser
from SummaryStats import SummaryStats

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
    v=float(totalAlt+1)/float(totalAlt+totalRef+2)
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

def runVariant(stan,ID,fields,numSamples,Lambda):
    # Write inputs file for STAN
    if(len(fields)<10): return (None,None)
    writeInputsFile(fields,INPUT_FILE)
    writeInitializationFile(fields,INIT_FILE)

    # Run STAN model
    if(DEBUG):
        cmd=stan.getCmd(WARMUP,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,
                        INIT_FILE)
        exit(cmd)
    stan.run(WARMUP,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE)

    # Parse output
    parser=StanParser(OUTPUT_TEMP)
    (thetaMedian,CI_left,CI_right)=parser.getMedianAndCI(1-ALPHA/2.0,"theta")
    p=parser.getVariable("p")
    q=parser.getVariable("q")
    (thetaMedian,thetaMean,thetaSD,thetaMin,thetaMax)=\
        parser.getSummary("theta")
    (pMedian,pMean,pSD,pMin,pMax)=parser.getSummary("p")
    (qMedian,qMean,qSD,qMin,qMax)=parser.getSummary("q")
    thetaVar=thetaSD*thetaSD; qVar=qSD*qSD; pVar=pSD*pSD
    Pleft=parser.getLeftTail("theta",1.0/Lambda)
    Pright=parser.getRightTail("theta",Lambda)
    Preg=max(Pleft,Pright)
    print(ID,round(thetaMedian,3),round(thetaVar,3),round(CI_left,3),
          round(CI_right,3),round(Preg,3),round(pMedian,3),round(pVar,5),
          round(qMedian,3),round(qVar,5),sep="\t")
    thetas=parser.getVariable("theta")
    thetas.sort(key=lambda x: x)
    return thetas

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:t:")
if(len(args)!=5):
    exit(ProgramName.get()+" [-s stanfile] [-t thetafile] <model> <min-effect> <input.txt> <#MCMC-samples> <firstVariant-lastVariant>\n   -s = save raw STAN file\n   -t = save theta samples\n   variant range is zero-based and inclusive\n   min-effect (lambda) must be >= 1\n")
(model,Lambda,inFile,numSamples,numVariants)=args
stanFile=None
thetaFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
    if(key=="-t"): thetaFile=value
if(not rex.find("(\d+)-(\d+)",numVariants)):
    exit(numVariants+": specify range of variants: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])
Lambda=float(Lambda)
if(Lambda<1): raise Exception("Min-effect must be >= 1")
THETA=None
if(thetaFile is not None): THETA=open(thetaFile,"wt")

stan=Stan(model)

# Process all input lines, each line = one variant (one MCMC run)
thetaIndex=None
variantIndex=0
print("ID\tTheta\tVar\tLeftCI\tRightCI\tPreg\tp\tvar(p)\tq\tvar(q)")

with open(inFile,"rt") as IN:
    for line in IN:
        # Check whether this variant is in the range to be processed
        if(variantIndex<firstIndex):
            variantIndex+=1
            continue
        elif(variantIndex>lastIndex): break
        fields=line.rstrip().split()
        ID=fields[0]
        thetas=runVariant(stan,ID,fields,numSamples,Lambda)
        variantIndex+=1
        if(thetas is None): continue
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

