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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <data.txt> <chrom-or-dot> <threshold>\n")
(infile,chrom,threshold)=sys.argv[1:]
threshold=int(threshold)

DEBUG=open("debug.txt","wt")
m={}
with open(infile,"rt") as IN:
    header=IN.readline()
    for line in IN:
        fields=line.rstrip().split("\t")
        n=len(fields)
        if(n!=14): continue
        (ID,v,output3_ref,output3_alt,output2_ref,output2_alt,output1_ref,
         output1_alt,input3_ref,input3_alt,input2_ref,input2_alt,input1_ref,
         input1_alt)=fields
        if(chrom!="."):
            if(not rex.find("(.*)@",ID)):
                raise Exception("Can't parse variant ID: "+ID)
            #if(rex[1]!=chrom): continue
        v=float(v)
        input1_alt=int(input1_alt); input1_ref=int(input1_ref)
        input2_alt=int(input2_alt); input2_ref=int(input2_ref)
        input3_alt=int(input3_alt); input3_ref=int(input3_ref)
        total=input3_ref+input3_alt+input2_ref+input2_alt+input1_ref+input1_alt
        if(total<threshold): continue
        alt=input1_alt+input2_alt+input3_alt
        ref=input1_ref+input2_ref+input2_ref
        p=float(alt)/float(alt+ref)
        if(m.get(v) is None): m[v]=[]
        m[v].append(p)
        zeros=0
        if(input1_alt==0 or input2_alt==0 or input3_alt==0 or
           input1_ref==0 or input2_ref==0 or input3_ref==0):
            zeros=1
        print(v,p,zeros,sep="\t")
        #if(v==0.5 and p<0.2):
        #if(v==0.5):
        if(True):
            print(ID,v,p,input1_alt,input1_ref,input2_alt,input2_ref,
                  input3_alt,input3_ref,sep="\t",file=DEBUG)

exit(0)
#========================= NOT RUN: ========================
V=m.keys()
V=[x for x in V]
V.sort(key=lambda x: x)
for v in V:
    array=m[v]
    (mean,SD,Min,Max)=SummaryStats.summaryStats(array)
    var=SD*SD
    if(v==0.5): print(v,var)    
    
