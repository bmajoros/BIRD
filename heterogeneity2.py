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

MIN=3

def printRecs(recs):
    for rec in recs:
        (x,mean,left,right)=rec
        label="="
        if(x<=left): label="<"
        elif(x>=right): label=">"
        print(label,x,mean,"("+str(left)+","+,str(right)+"),sep="\t")
    print()

def process(recs):
    n=len(recs)
    if(n<MIN): return
    numLeft=0; numRight=0
    for rec in recs:
        (x,left,right)=rec
        if(x<=left): numLeft+=1
        elif(x>=right): numRight+=1
    if(numLeft>0 and numRight>0):
        printRecs(recs)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <output-of-heterogeneity.py.txt>\n")
(infile,)=sys.argv[1:]

prevID=None
recs=[]
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=6): continue
        (ID,inout,theoretical,mean,left,right)=fields
        if(ID!=prevID):
            if(len(recs)>0):
                process(recs)
                recs=[]
        if(inout=="OUT"):
            recs.append([theoretical,mean,left,right])
        prevID=ID

