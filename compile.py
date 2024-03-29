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
import os
import ProgramName
from Rex import Rex
rex=Rex()

HERE="/hpc/group/majoroslab/BIRD/git/"
STAN="/hpc/group/majoroslab/cmdstan-2.27.0"

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <model>\n")
(model,)=sys.argv[1:]

if(rex.find(".stan$",model)):
    exit("do not include file extension in model name")
os.system("cd "+STAN+ "; rm -f "+HERE+model+".o; make "+HERE+model)
