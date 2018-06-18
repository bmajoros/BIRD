#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
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
from Pipe import Pipe
from SummaryStats import SummaryStats
import getopt
import SumLogProbs

WANT_BAYES_FACTOR=False
BAYES_NULL="prior-1"
BAYES_LESS="prior-less"
BAYES_GREATER="prior-greater"
WARMUP=1000
ALPHA=0.05
STDERR=TempFilename.generate(".stantimes")
INPUT_FILE=TempFilename.generate(".staninputs")
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

def getThetaIndex(fields):
    numFields=len(fields)
    thetaIndex=-1
    pIndex=-1
    qIndex=-1
    for i in range(7,numFields):
        if(fields[i]=="theta"): thetaIndex=i
        elif(fields[i]=="p"): pIndex=i
        elif(fields[i]=="q"): qIndex=i
    return (thetaIndex,pIndex,qIndex)

def writeToFile(fields,OUT):
    numFields=len(fields)
    for i in range(7,numFields):
        print(fields[i],end="",file=OUT)
        if(i<numFields-1): print("\t",end="",file=OUT)
    print(file=OUT)

def writeReadCounts(fields,start,numReps,varName,OUT,packArrays):
    if(packArrays):
        print(varName,"<- array(c(",file=OUT,end="")
    else:
        print(varName,"<- c(",file=OUT,end="")
    for rep in range(numReps):
        print(fields[start+rep*2],file=OUT,end="")
        if(rep+1<numReps): print(",",file=OUT,end="")
    if(packArrays):
        print("),dim="+str(numReps)+")",file=OUT)
    else:
        print(")",file=OUT)

def writeInputsFile(fields,filename,swap,packArrays):
    v=fields[0]; DNAreps=int(fields[4])
    rnaIndex=5+2*DNAreps
    RNAreps=int(fields[rnaIndex])
    OUT=open(filename,"wt")
    print("v <-",v,file=OUT)
    print("N_DNA <-",str(DNAreps),file=OUT)
    if(not swap):
        writeReadCounts(fields,6,DNAreps,"a",OUT,packArrays) # alt
        writeReadCounts(fields,5,DNAreps,"b",OUT,packArrays) # ref
    else:
        writeReadCounts(fields,5,DNAreps,"a",OUT,packArrays) # alt
        writeReadCounts(fields,6,DNAreps,"b",OUT,packArrays) # ref
    print("N_RNA <-",str(RNAreps),file=OUT)
    if(not swap):
        writeReadCounts(fields,rnaIndex+2,RNAreps,"k",OUT,packArrays) # alt
        writeReadCounts(fields,rnaIndex+1,RNAreps,"m",OUT,packArrays) # ref
    else:
        writeReadCounts(fields,rnaIndex+1,RNAreps,"k",OUT,packArrays) # alt
        writeReadCounts(fields,rnaIndex+2,RNAreps,"m",OUT,packArrays) # ref
    OUT.close()

def getMeanSD(thetas):
    (mean,SD,min,max)=SummaryStats.summaryStats(thetas)
    return (mean,SD)

def getMedian(thetas):
    # Precondition: thetas is already sorted
    return thetas[int(len(thetas)/2)]

def getCredibleInterval(thetas,alpha):
    halfAlpha=alpha/2.0
    n=len(thetas)
    leftIndex=int(halfAlpha*n)
    rightIndex=n-leftIndex
    left=thetas[leftIndex+1]
    right=thetas[rightIndex-1]
    return (left,right)

def isSignificant(thetas):
    n=len(thetas)
    median=getMedian(thetas)
    oneOrLess=0
    oneOrGreater=0
    for i in range(n):
        if(thetas[i]<=1): oneOrLess+=1
        if(thetas[i]>=1): oneOrGreater+=1
    leftP=float(oneOrLess)/float(n)
    rightP=float(oneOrGreater)/float(n)
    P=leftP if median>1 else rightP
    return P<ALPHA

def computeThetaDirectly(fields,pIndex,qIndex):
    if(pIndex<0 or qIndex<0): return 1 # beta-binomial model
    p=float(fields[pIndex])
    q=float(fields[qIndex])
    theta=q/p/((1-q)/(1-p))
    return theta

def runVariant(model,fields,thin,numSamples,outfile,swap):
    # Write inputs file for STAN
    if(len(fields)<10): return (None,None)
    writeInputsFile(fields,INPUT_FILE,swap,False)

    # Run STAN model
    cmd=model+" sample thin="+str(thin)+\
        " num_samples="+numSamples+\
        " num_warmup="+str(WARMUP)+\
        " data file="+INPUT_FILE+\
        " output file="+OUTPUT_TEMP+" > "+STDERR
    #print(cmd)
    os.system(cmd)

    # Parse MCMC output
    c1Index=None; c2Index=None; c3Index=None
    thetas=[]
    OUT=open(outfile,"wt")
    with open(OUTPUT_TEMP,"rt") as IN:
        for line in IN:
            if(len(line)==0 or line[0]=="#"): continue
            fields=line.rstrip().split(",")
            numFields=len(fields)
            if(numFields>0 and fields[0]=="lp__"):
                printFields(fields,OUT)
                (thetaIndex,pIndex,qIndex)=getThetaIndex(fields)
                c1Index=getFieldIndex("c1",fields)
                c2Index=getFieldIndex("c2",fields)
                c3Index=getFieldIndex("c3",fields)
            else:
                writeToFile(fields,OUT)
                theta=None
                if(thetaIndex>=0):
                    theta=float(fields[thetaIndex])
                else:
                    theta=computeThetaDirectly(fields,pIndex,qIndex)
                thetas.append(theta)
    OUT.close()
    autoCorLen=-1
    thetas.sort(key=lambda x: x)
    return(thetas,autoCorLen)

def getStdError(thetas):
    logThetas=[]
    for theta in thetas:
        logThetas.append(math.log(theta)/math.log(2))
    (mean,SD,min,max)=SummaryStats.summaryStats(logThetas)
    return SD/math.sqrt(len(thetas))

def summarize(thetas,autoCorLen,fields):
    n=len(thetas)
    mode=0
    (mean,sd)=getMeanSD(thetas)
    #stdError=getStdError(thetas) #sd/math.sqrt(n)
    median=getMedian(thetas)
    (CI_left,CI_right)=getCredibleInterval(thetas,ALPHA)
    oneOrLess=0
    oneOrGreater=0
    for i in range(n):
        if(thetas[i]<=1): oneOrLess+=1
        if(thetas[i]>=1): oneOrGreater+=1
    leftP=float(oneOrLess)/float(n)
    rightP=float(oneOrGreater)/float(n)
    P=leftP if leftP<rightP else rightP
    BF_less=1
    BF_greater=1
    if(WANT_BAYES_FACTOR):
        if(USE_BRIDGE_SAMPLING):
            writeInputsFile(fields,INPUT_FILE,False,True)
            (BF_less,BF_greater)=bridgeSampler(INPUT_FILE)
        else:
            BF_less=computeBayesFactor(BAYES_LESS,fields)
            BF_greater=computeBayesFactor(BAYES_GREATER,fields)
    print(median,CI_left,CI_right,BF_less,BF_greater,sep="\t")

def computeBayesFactor(model,field):
    nullValue=runBayesFactorModel(BAYES_NULL,fields)
    if(nullValue is None): return 1
    altValue=runBayesFactorModel(model,fields)
    if(altValue is None): return 1
    BF=math.exp(altValue-nullValue)
    return BF

def runBayesFactorModel(model,fields):
    thin=1
    numSamples=1000

    # Write inputs file for STAN
    if(len(fields)<10): return None
    writeInputsFile(fields,INPUT_FILE,False,False)

    # Run STAN model
    cmd=model+" sample thin="+str(thin)+\
        " num_samples="+str(numSamples)+\
        " num_warmup="+str(WARMUP)+\
        " data file="+INPUT_FILE+\
        " output file="+OUTPUT_TEMP+" > "+STDERR
    #print(cmd)
    os.system(cmd)

    # Parse MCMC output
    logLikelihoods=[]
    LLindex=None
    with open(OUTPUT_TEMP,"rt") as IN:
        for line in IN:
            if(len(line)==0 or line[0]=="#"): continue
            fields=line.rstrip().split(",")
            numFields=len(fields)
            if(numFields>0 and fields[0]=="lp__"):
                LLindex=getFieldIndex("LL",fields)
            else:
                if(LLindex is None): raise Exception("No LL field in output")
                else:
                    LL=float(fields[LLindex])
                    logLikelihoods.append(LL)
    n=len(logLikelihoods)
    if(n<1):
        raise Exception("No likelihoods found in STAN output for model "+model)
    logLikelihood=SumLogProbs.sumLogProbs(logLikelihoods)-math.log(n)
    return logLikelihood

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:")
if(len(args)!=5):
    exit(ProgramName.get()+" [-s file] <model> <input.txt> <output.txt> <#MCMC-samples> <firstVariant-lastVariant>\n   -s = save raw STAN file\n   variant range is zero-based and inclusive\n")
(model,inputFile,outfile,numSamples,numVariants)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
if(not rex.find("(\d+)-(\d+)",numVariants)):
    exit(numVariants+": specify range of variants: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])

# Process all input lines, each line = one variant (one MCMC run)
thetaIndex=None; pIndex=None; qIndex=None
variantIndex=0
with open(inputFile,"rt") as IN:
    for line in IN:
        # Check whether this variant is in the range to be processed
        if(variantIndex<firstIndex):
            variantIndex+=1
            continue
        elif(variantIndex>lastIndex): break
        fields=line.rstrip().split()
        (thetas,autoCorLen)=runVariant(model,fields,1,numSamples,outfile,
                                       False)
        if(thetas is None): continue
        summarize(thetas,autoCorLen,fields)
        variantIndex+=1
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)

