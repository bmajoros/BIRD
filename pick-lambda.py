#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2018 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import os

BIRD=os.environ['BIRD']
if(BIRD is None): exit("Please set the $BIRD environment variable")
TABLE=BIRD+"/lambda.txt"
FDRS=(0.1, 0.05, 0.01)
DEPTHS=(50000, 10000, 1000, 100)
THETAS=(2, 1.5, 1.25, 1/1.25, round(1/1.5,2), 1/2)
MAFS=(0.5, 0.10, 0.05, 0.01)

class Record:
    def __init__(this,Lambda,fdr,fdrI,theta,thetaI,maf,mafI,depth,depthI):
        this.Lambda=Lambda
        this.fdr=fdr
        this.fdrI=fdrI
        this.theta=theta
        this.thetaI=thetaI
        this.maf=maf
        this.mafI=mafI
        this.depth=depth
        this.depthI=depthI
    def manhattan(this,fdrI,thetaI,mafI,depthI):
        return abs(this.fdrI-fdrI)+abs(this.thetaI-thetaI)+\
            abs(this.mafI-mafI)+abs(this.depthI-depthI)
    def matches(this,fdr,theta,maf,depth):
        #print(this.fdr,fdr,"\t",this.theta,theta,"\t",this.maf,maf,"\t",
        #      this.depth,depth)
        return this.fdr==fdr and this.theta==theta and this.maf==maf and\
            this.depth==depth

def loadTable(filename):
    records=[]
    with open(filename,"rt") as IN:
        IN.readline() # header
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)==0): continue
            if(len(fields)!=7): raise Exception("can't parse table")
            (Lambda,alpha,fdr,power,theta,maf,depth)=fields
            Lambda=float(Lambda); alpha=float(alpha); theta=float(theta)
            maf=float(maf); depth=int(depth)
            fdrI=findNearest(alpha,FDRS)[1]
            thetaI=findNearest(theta,THETAS)[1]
            mafI=findNearest(maf,MAFS)[1]
            depthI=findNearest(depth,DEPTHS)[1]
            rec=Record(Lambda,alpha,fdrI,theta,thetaI,
                       maf,mafI,depth,depthI)
            records.append(rec)
    return records

def findNearest(value,legal):
    n=len(legal)
    best=None
    bestDiff=None
    bestI=None
    for i in range(n):
        x=legal[i]
        diff=abs(value-x)
        if(best is None or diff<bestDiff):
            best=x
            bestDiff=diff
            bestI=i
    return (best,bestI)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit("\n"+ProgramName.get()+" <FDR> <effect> <maf> <depth>\n\
  FDR is one of 0.1, 0.05, 0.1\n\
  effect is ratio of transcriptional reates for alt versus ref allele\n\
")
(fdr,theta,maf,depth)=sys.argv[1:]
fdr=float(fdr)
theta=float(theta)
maf=float(maf)
depth=int(depth)

# Check that parms are legal
if(fdr!=0.1 and fdr!=0.05 and fdr!=0.01):
    exit("FDR must be one of 0.1, 0.05, 0.01")
if(theta<0): exit("effect must be > 0")
#if(theta>1): theta=1.0/theta
if(maf<=0 or maf>=1): exit("maf must be between 0 and 1")
if(depth<=0): exit("depth must be > 0")

# Find closest simulation parameters
fdrI=findNearest(fdr,FDRS)[1]
(theta,thetaI)=findNearest(theta,THETAS)
(maf,mafI)=findNearest(maf,MAFS)
(depth,depthI)=findNearest(depth,DEPTHS)
#print("fdr=",fdr,"theta=",theta,"maf=",maf,"depth=",depth)

table=loadTable(TABLE)
minDist=None
best=None
bestDist=None
for rec in table:
    dist=rec.manhattan(fdrI,thetaI,mafI,depthI)
    if(minDist is None or dist<minDist):
        best=rec
        bestDist=dist
print(best.Lambda)


