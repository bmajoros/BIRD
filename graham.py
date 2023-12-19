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
import ProgramName
from Pipe import Pipe
from Rex import Rex
rex=Rex()
from SummaryStats import SummaryStats

def loadBIRD(filename):
    byVariant={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)<5): continue
            (ID,theta,left,right,Preg)=fields
            byVariant[ID]=float(Preg)
    return byVariant

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <inputs.txt> <outputs.txt> <threshold>\n")
(inputsFile,outputsFile,threshold)=sys.argv[1:]
threshold=int(threshold)

BIRD=loadBIRD(outputsFile)
DEBUG=open("debug.txt","wt")
m={}
with open(inputsFile,"rt") as IN:
    header=IN.readline()
    for line in IN:
        fields=line.rstrip().split("\t")
        n=len(fields)
        if(n!=21): continue
        ID=fields[0]
        dnaRef=[fields[2],fields[4],fields[6],fields[8],fields[10]]
        dnaAlt=[fields[3],fields[5],fields[7],fields[9],fields[11]]
        rnaRef=[fields[13],fields[15],fields[17],fields[19]]
        rnaAlt=[fields[14],fields[16],fields[18],fields[20]]
        dnaAltSum=0; dnaRefSum=0; zeros=0
        for i in range(5):
            dnaRef[i]=int(dnaRef[i]); dnaAlt[i]=int(dnaAlt[i])
            dnaRefSum+=dnaRef[i];     dnaAltSum+=dnaAlt[i]
            if(dnaRef[i]==0 or dnaAlt[i]==0): zeros=1
        if(dnaAltSum==0 or dnaRefSum==0): continue
        dnaSum=dnaAltSum+dnaRefSum
        if(dnaSum<threshold): continue
        p=float(dnaAltSum)/float(dnaAltSum+dnaRefSum)
        v=0.5
        if(m.get(v) is None): m[v]=[]
        m[v].append(p)
        print(v,p,zeros,sep="\t")
        #if(v==0.5 and p<0.2):
        #if(v==0.5):
        significant="no"
        if(ID in BIRD and BIRD[ID]>0.9): significant="yes"
        if(True):
            print(ID,v,p,dnaAlt[0],dnaRef[0],dnaAlt[1],dnaRef[1],
                  dnaAlt[2],dnaRef[2],dnaAlt[3],dnaRef[3],
                  dnaAlt[4],dnaRef[4],
                  significant,sep="\t",file=DEBUG)

    
