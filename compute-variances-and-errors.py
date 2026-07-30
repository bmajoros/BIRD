#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import statistics
import math
import ProgramName
from Rex import Rex
rex=Rex()

class Variant:
    def __init__(self,ID,groups,p_bar,q_bar):
        self.ID=ID
        self.groups=groups
        self.p_bar=p_bar # average AF across groups (pools or reps) in DNA
        self.q_bar=q_bar # ditto, but for RNA instead of DNA

class Group:
    def __init__(self,dna,rna):
        self.dna=dna # [ref,alt]
        self.rna=rna # [ref,alt]

def readData(filename):
    variants=[]
    with open(filename,"rt") as IN:
        header=IN.readline()
        rex.findOrDie("variants=(\d+) groups=(\d+) heterogeneity=(\S+) theta=(\S+)",header)
        numVariants=int(rex[1])
        numGroups=int(rex[2])
        heterogeneity=float(rex[3])
        theta=float(rex[4])
        for line in IN:
            v=parseVariant(line,numGroups)
            variants.append(v)
    return (variants,numVariants,numGroups,heterogeneity,theta)

def parseVariant(line,numGroups):
    fields=line.rstrip().split()
    ID=fields[0]
    rex.findOrDie("p_bar=(\S+)",fields[1])
    p_bar=float(rex[1])
    rex.findOrDie("q_bar=(\S+)",fields[2])
    q_bar=float(rex[1])
    nextField=3
    groups=[]
    for i in range(numGroups):
        dna=parseCounts(fields[nextField])
        nextField+=1
        rna=parseCounts(fields[nextField])
        nextField+=1
        group=Group(dna,rna)
        groups.append(group)
    return Variant(ID,groups,p_bar,q_bar)

def parseCounts(text):
    rex.findOrDie("na=(\d+),(\d+),(\S+)",text)
    ref=int(rex[1])
    alt=int(rex[2])
    pq=float(rex[3])
    return [ref,alt,pq]

def computeVarDNA(variants):
    values=[]
    for variant in variants:
        A=sum([group.dna[1] for group in variant.groups]) # sum alt counts
        values.append(A)
    return statistics.variance(values)

def computeVarRNA(variants):
    values=[]
    for variant in variants:
        A=sum([group.rna[1] for group in variant.groups]) # sum alt counts
        values.append(A)
    return statistics.variance(values)

def computeVarP(variants):
    values=[]
    squaredErrors=[]
    for variant in variants:
        R=sum([group.dna[0] for group in variant.groups]) # sum ref counts
        A=sum([group.dna[1] for group in variant.groups]) # sum alt counts
        p_hat=float(A)/float(A+R)
        variant.p_hat=p_hat
        values.append(p_hat)
        error=p_hat=variant.p_bar
        squaredErrors.append(error*error)
    var=statistics.variance(values)
    rmse=math.sqrt(sum(squaredErrors)/len(squaredErrors))
    return (var,rmse)

def computeVarQ(variants):
    values=[]
    squaredErrors=[]
    for variant in variants:
        R=sum([group.rna[0] for group in variant.groups]) # sum ref counts
        A=sum([group.rna[1] for group in variant.groups]) # sum alt counts
        q_hat=float(A)/float(A+R)
        variant.q_hat=q_hat
        values.append(q_hat)
        error=q_hat=variant.q_bar
        squaredErrors.append(error*error)
    var=statistics.variance(values)
    rmse=math.sqrt(sum(squaredErrors)/len(squaredErrors))
    return (var,rmse)

def computeVarTheta(variants,theta):
    values=[]
    squaredErrors=[]
    for variant in variants:
        p=variant.p_hat
        q=variant.q_hat
        theta_hat=q/(1-q)/(p/(1-p))
        values.append(theta_hat)
        error=theta_hat-theta
        squaredErrors.append(error*error)
    var=statistics.variance(values)
    rmse=math.sqrt(sum(squaredErrors)/len(squaredErrors))
    return (var,rmse)
    

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <counts.txt>\n")
(infile,)=sys.argv[1:]

# Parse filename
rex.findOrDie("pools(\d+)-theta(\S+)-reads(\d+)-noise(\S+).txt",infile)
numPools=rex[1]
theta=rex[2]
reads=rex[3]
noise=rex[4]

# Load variants
variants,numVariants,numGroups,heterogeneity,theta=readData(infile)

# Compute variance in DNA counts
varDNA=computeVarDNA(variants)
varRNA=computeVarRNA(variants)

# Compute variance and RMSE in estimates of p and q
(varP,rmseP)=computeVarP(variants)
(varQ,rmseQ)=computeVarQ(variants)

# Compute variance and RMSE in estimate of theta
(varTheta,rmseTheta)=computeVarTheta(variants,theta)

# Report stats
print("pools\ttheta\treads\tnoise\tvarDNA\tvarRNA\tvarP\tvarQ\trmseP\t"+\
      "rmseQ\tvarTheta\trmseTheta")
print(numPools,theta,reads,noise,round(varDNA,4),round(varRNA,4),
      round(varP,7),round(varQ,7),round(rmseP,5),round(rmseQ,5),
      round(varTheta,5),round(rmseTheta,5),sep="\t")

#print("pools=",numPools," theta=",theta," reads=",reads," noise=",noise,
#      " varDNA=",round(varDNA,2)," varRNA=",round(varRNA,3),
#      " varP=",round(varP,7)," varQ=",round(varQ,7),
#      " rmseP=",round(rmseP,4)," rmseQ=",round(rmseQ,4),
#      " varTheta=",round(varTheta,3)," rmseTheta=",round(rmseTheta,3),
#      sep="")

