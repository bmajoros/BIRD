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
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <debug.txt>\n")
(infile,)=sys.argv[1:]

prevPos=None; prevChrom=None
with open(infile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)!=10): continue
        (ID,v,estimate,alt1,ref1,alt2,ref2,alt3,ref3,significant)=fields
        v=float(v); estimate=float(estimate)
        #diff=(v-estimate)/v
        diff=estimate-v
        #if(diff<-1): print("XXX",v,estimate,sep="\t")
        if(not rex.find("(.*)@(\d+)",ID)):
            raise Exception("Can't parse: "+ID)
        chrom=rex[1]; pos=int(rex[2])
        if(chrom!=prevChrom): prevPos=None
        if(prevPos is not None):
            distance=pos-prevPos
            if(distance<0): raise Exception("not sorted")
            print(distance,diff,significant,sep="\t")
        prevPos=pos
        prevChrom=chrom

        


