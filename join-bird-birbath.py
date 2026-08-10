#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import ProgramName

def loadBird(poolsDir):
    pass

def processBirdBath(filename,BIRD):
    pass

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <birdbath.txt> <bird-pools-dir>\n")
(birdbathFilename,poolsDir)=sys.argv[1:]

BIRD=loadBird(poolsDir)
processBirdBath(birdbathFilename,BIRD)


