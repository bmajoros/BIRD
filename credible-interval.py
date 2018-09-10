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
import numpy as np

def toFloat(fields):
    array=[]
    for field in fields: array.append(float(field))
    return array

def analyze(fields,lower,upper):
    array=toFloat(fields)
    npArray=np.array(array)
    a=np.percentile(npArray,lower*100)
    b=np.percentile(npArray,upper*100)
    return (a,b)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <thetas.txt> <alpha>\n")
(filename,alpha)=sys.argv[1:]
alpha=float(alpha)
lower=alpha/2.0
upper=1.0-lower

with open(filename,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        (a,b)=analyze(fields,lower,upper)
        print(a,b,sep="\t")


