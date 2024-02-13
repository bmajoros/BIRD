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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <in:pooled.essex> <in:collapsed.essex> <threshold> <out:pooled.essex> <out:collapsed.essex>\n")
(pooledInfile,collapsedInfile,threshold,pooledOutfile,collapsedOutfile)=\
    sys.argv[1:]
threshold=float(threshold)

OUT_POOLED=open(pooledOutfile,"wt")
OUT_COLLAPSED=open(collapsedOutfile,"wt")
pooledParser=PooledParser(pooledInfile)
collapsedParser=PooledParser(collapsedInfile)
while(True):
    varP=pooledParser.nextVariant()
    varC=collapsedParser.nextVariant()
    if(varP is None or varC is None): break
    if(varP.ID!=varC.ID): raise Exception("Files not in sync")
    freqs=varC.getFreqs()
    if(len(freqs)!=1): raise Exception("Freqs array has unexpected size")
    maf=freqs[0]
    #print(maf,threshold,maf<=threshold,sep="\t")
    if(maf>threshold): continue
    textP=varP.print()
    textC=varC.print()
    print(textP,file=OUT_POOLED,flush=True)
    print(textC,file=OUT_COLLAPSED,flush=True)
OUT_POOLED.close(); OUT_COLLAPSED.close()






