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
from PooledParser import PooledParser
from PooledVariant import PooledVariant

def loadTruth(infile):
    truth={}
    with open(infile,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): continue
            (ID,theta)=fields
            theta=float(theta)
            truth[ID]=theta
    return truth

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <in.essex> <truth.txt> <dna_freq> <out.essex>\n")
(infile,truthFile,p,outfile)=sys.argv[1:]
p=float(p) # DNA freq

truth=loadTruth(truthFile)
OUT=open(outfile,"wt")
parser=PooledParser(infile)
while(True):
    var=parser.nextVariant()
    if(var is None): break
    theta=truth[var.ID] # effect size
    q=theta*p/(1-p+theta*p) # RNA freq
    var.forceEqualFreqs(p,q)
    text=var.print()
    print(text,file=OUT,flush=True)
OUT.close()




