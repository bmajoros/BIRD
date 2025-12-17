#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import math
import ProgramName

MIN=10

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <pool-genotypes.txt> <bird-input.txt>\n")
(genotypeFile,countsFile)=sys.argv[1:]

freqs=dict()
with open(genotypeFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=6): continue
        (chrom,pos,label,ref,alt,c)=fields
        ID=chrom+":"+pos+":"+ref+":"+alt
        c=int(c)
        if(c==0 or c==10): continue
        freq=float(c)/10.0
        freqs[ID]=freq

numIn=0; numOut=0
with open(countsFile,"rt") as IN:
    for line in IN:
        if("NA" in line): continue
        fields=line.rstrip().split()
        if(len(fields)!=15): continue
        DNA_START=2
        RNA_START=9
        START=RNA_START
        DNA=([fields[START],fields[START+1]],
             [fields[START+2],fields[START+3]],
             [fields[START+4],fields[START+5]])
        if(int(DNA[0][0])+int(DNA[0][1])<MIN): continue
        if(int(DNA[1][0])+int(DNA[1][1])<MIN): continue
        if(int(DNA[2][0])+int(DNA[2][1])<MIN): continue
        ID=fields[0]
        freq=freqs.get(ID,None)
        if(freq is None): continue
        for dna in DNA:
            (ref,alt)=(int(dna[0]),int(dna[1]))
            n=float(alt+ref)
            p=float(alt)/n
            SE=math.sqrt(p*(1-p)/n)
            if(SE==0.0): continue
            interval=(p-1.96*SE,p+1.96*SE)
            label="IN" if (freq>=interval[0] and freq<=interval[1]) else "OUT"
            print(ID,label,freq,round(p,4),round(interval[0],2),round(interval[1],2),sep="\t")
            if(label=="IN"): numIn+=1
            else: numOut+=1
print(numIn,"IN",float(numIn)/float(numIn+numOut))
print(numOut,"OUT",float(numOut)/float(numIn+numOut))

