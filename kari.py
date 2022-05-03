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

def getNaiveEstimate(ref1,alt1,ref2,alt2,ref3,alt3):
    ref1=float(ref1); alt1=float(alt1)
    ref2=float(ref2); alt2=float(alt2)
    ref3=float(ref3); alt3=float(alt3)
    ref=ref1+ref2+ref3
    alt=alt1+alt2+alt3
    af=alt/(alt+ref)
    return round(af,3)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <data.txt>\n")
(infile,)=sys.argv[1:]

with open(infile,"rt") as IN:
    header=IN.readline()
    for line in IN:
        fields=line.rstrip().split("\t")
        n=len(fields)
        if(n!=14): continue
        (ID,v,output3_ref,output3_alt,output2_ref,output2_alt,output1_ref,
         output1_alt,input3_ref,input3_alt,input2_ref,input2_alt,input1_ref,
         input1_alt)=fields
        fields=(ID,3,input3_ref,input3_alt,input2_ref,input2_alt,input1_ref,
                input1_alt,3,output3_ref,output3_alt,output2_ref,output2_alt,
                output1_ref,output1_alt)
        naiveQ=getNaiveEstimate(output3_ref,output3_alt,output2_ref,
                                output2_alt,output1_ref,output1_alt)
        naiveP=getNaiveEstimate(input3_ref,input3_alt,input2_ref,
                                input2_alt,input1_ref,input1_alt)
        fields=[str(x) for x in fields]
        line="\t".join(fields)
        OUT=open("temp.txt","wt")
        print(line,file=OUT)
        OUT.close()
        cmd="git/bird2.py git/BIRD 1.2 temp.txt 1000 0-1"
        pipe=Pipe(cmd)
        header=pipe.readline()
        output=None
        while(True):
            line=pipe.readline()
            if(line is None): break
            if(rex.find("Informational",line)): continue
            if(rex.find("Exception",line)): continue
            if(rex.find("warning",line)): continue
            output=line
        pipe.close()
        print(v,naiveP,naiveQ,output,sep="\t")

