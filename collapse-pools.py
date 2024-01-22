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
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <in.essex> <out.essex>\n")
(infile,outfile)=sys.argv[1:]

OUT=open(outfile,"wt")
parser=PooledParser(infile)
while(True):
    var=parser.nextVariant()
    if(var is None): break
    var=var.collapse()
    text=var.print()
    print(text,file=OUT,flush=True)
OUT.close()




