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

def writeInitializationFile(stan,filename):
    OUT=open(filename,"wt")
    print("theta <- 1",file=OUT)
    print("p <- 0.5",file=OUT)
    OUT.close()

def writeInputsFile(stan,dna_ref,dna_alt,rna_ref,rna_alt,b1,b2,filename):
    OUT=open(filename,"wt")
    print("dna_ref <-",dna_ref,file=OUT)
    print("dna_alt <-",dna_alt,file=OUT)
    print("rna_ref <-",rna_ref,file=OUT)
    print("rna_alt <-",rna_alt,file=OUT)
    print("beta1 <-",b1,file=OUT)
    print("beta2 <-",b2,file=OUT)
    OUT.close()

def runVariant(stan,dna_ref,dna_alt,rna_ref,rna_alt,b1,b2,numSamples):
    # Write inputs file for STAN
    writeInputsFile(stan,dna_ref,dna_alt,rna_ref,rna_alt,b1,b2,INPUT_FILE)
    writeInitializationFile(stan,INIT_FILE)

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
if(len(args)!=9):
    exit(ProgramName.get()+" [-s stanfile] [-t thetafile] <model> <min-effect> <#MCMC-samples> <dna_ref> <dna_alt> <rna_ref> <rna_alt> <beta-parm1> <beta-parm2>\n   -s = save raw STAN file\n   -t = save theta samples\n   min-effect (lambda) must be >= 1\n")
(model,minEffect,numSamples,dna_ref,dna_alt,rna_ref,rna_alt,b1,b2)=args
stanFile=None
thetaFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
    if(key=="-t"): thetaFile=value
minEffect=float(minEffect)
dna_ref=int(dna_ref); dna_alt=int(dna_alt)
rna_ref=int(rna_ref); rna_alt=int(rna_alt)
b1=float(b1)
b2=float(b2)
if(minEffect<1): raise Exception("Min-effect must be >= 1")
THETA=None
if(thetaFile is not None): THETA=open(thetaFile,"wt")
stan=Stan(model)

(thetas,stanParser)=runVariant(stan,dna_ref,dna_alt,rna_ref,rna_alt,
                               b1,b2,numSamples)
#summarize(stanParser,thetas,variant.ID,minEffect)
if(THETA is not None):
    for i in range(len(thetas)):
        print(thetas[i],file=THETA)
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)
if(THETA is not None): THETA.close()

