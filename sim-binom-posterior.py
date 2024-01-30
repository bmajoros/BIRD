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
import numpy as np

SAMPLES=10000

def sampleBeta(k,m,mode,conc):
    alpha=mode*(conc-2)+1
    beta=(1-mode)*(conc-2)+1
    s=np.random.beta(k+alpha,m+beta)
    return s

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" k1 n1 k2 n2 w1 w2 c\n")
(k1,n1,k2,n2,w1,w2,c)=sys.argv[1:]
k1=int(k1); k2=int(k2); n1=int(n1); n2=int(n2)
w1=float(w1); w2=float(w2); c=float(c)

m1=n1-k1
m2=n2-k2
for i in range(SAMPLES):
    # Model with pools:
    q1=sampleBeta(k1,m1,w1,c)
    q2=sampleBeta(k2,m2,w2,c)
    ave1=(q1+q2)/2

    # Collapsed model:
    ave2=sampleBeta(k1+k2,m1+m2,(w1+w2)/2,c)
    print(q1,q2,ave1,ave2,sep="\t")
    



