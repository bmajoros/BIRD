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
import ProgramName

class Record:
    def __init__(this,id,theta,left,right,Preg):
        this.id=id
        this.theta=theta
        this.left=left
        this.right=right
        this.Preg=Preg
        this.fdr=None
        this.sig="nonsignificant"
    def print(this):
        print(this.id,this.theta,this.left,this.right,this.Preg,this.fdr,
              this.sig,sep="\t")

def loadPredictions(filename):
    recs=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=5): raise Exception("wrong number of fields")
            (variantID,theta,left,right,Preg)=fields
            Preg=float(Preg)
            rec=Record(variantID,theta,left,right,Preg)
            recs.append(rec)
    return recs

def justin(recs,ALPHA):
    n=len(recs)
    sumP=0
    keep=0
    for i in range(n):
        sumP+=1-recs[i].Preg
        recs[i].fdr=sumP/float(i+1)
        if(recs[i].fdr>ALPHA): break
        else: recs[i].sig="SIGNIFICANT"
        keep=i
    return keep

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <predictions.txt> <fdr>\n")
(infile,fdr)=sys.argv[1:]
fdr=float(fdr)

preds=loadPredictions(infile)
preds.sort(key=lambda x: 1-x.Preg)
keep=justin(preds,fdr)
for rec in preds: rec.print()
    
