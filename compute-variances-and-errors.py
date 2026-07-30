#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import statistics
import ProgramName
from Rex import Rex
rex=Rex()

class Variant:
    def __init__(self,ID,groups):
        self.ID=ID
        self.groups=groups

class Group:
    def __init__(self,dna,rna):
        self.dna=dna # [ref,alt]
        self.rna=rna # [ref,alt]

def readData(filename):
    variants=[]
    with open(filename,"rt") as IN:
        header=IN.readline()
        rex.findOrDie("variants=(\d+) groups=(\d+) heterogeneity=(\S+)",header)
        numVariants=int(rex[1])
        numGroups=int(rex[2])
        heterogeneity=float(rex[3])
        for line in IN:
            v=parseVariant(line,numGroups)
            variants.append(v)
    return (variants,numVariants,numGroups,heterogeneity)

def parseVariant(line,numGroups):
    fields=line.rstrip().split()
    ID=fields[0]
    nextField=1
    groups=[]
    for i in range(numGroups):
        dna=parseCounts(fields[nextField])
        nextField+=1
        rna=parseCounts(fields[nextField])
        nextField+=1
        group=Group(dna,rna)
        groups.append(group)
    return Variant(ID,groups)

def parseCounts(text):
    rex.findOrDie("na=(\d+),(\d+)",text)
    ref=int(rex[1])
    alt=int(rex[2])
    return [ref,alt]

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
    for variant in variants:
        R=sum([group.dna[0] for group in variant.groups]) # sum ref counts
        A=sum([group.dna[1] for group in variant.groups]) # sum alt counts
        p_hat=float(A)/float(A+R)
        values.append(p_hat)
    return statistics.variance(values)

def computeVarQ(variants):
    values=[]
    for variant in variants:
        R=sum([group.rna[0] for group in variant.groups]) # sum ref counts
        A=sum([group.rna[1] for group in variant.groups]) # sum alt counts
        p_hat=float(A)/float(A+R)
        values.append(p_hat)
    return statistics.variance(values)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <counts.txt>\n")
(infile,)=sys.argv[1:]

# Load variants
variants,numVariants,numGroups,heterogeneity=readData(infile)

# Compute variance in DNA counts
varDNA=computeVarDNA(variants)
varRNA=computeVarRNA(variants)
print(varDNA,varRNA)

# Compute variance in estimates of p and q
varP=computeVarP(variants)
varQ=computeVarQ(variants)
print(varP,varQ)
