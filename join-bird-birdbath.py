#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
from pathlib import Path
import ProgramName
from Rex import Rex

rex=Rex()

class Variant:
    def __init__(self,ID):
        self.ID=ID
        self.pools=[]
        
def loadBird(poolsDir):
    variants=dict()
    files=[str(f) for f in Path(poolsDir).iterdir() if f.is_file()]
    for filename in files:
        if(not rex.find("out(\d+).txt",filename)): continue
        readBird(filename,variants,rex[1])
    return variants

def readBird(filename,variants,poolID):
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=5): continue
            (ID,theta,left,roght,Preg)=fields
            if(not rex.find("^chr",ID)): continue
            variant=variants.get(ID,None)
            if(variant is None):
                variant=variants[ID]=Variant(ID)
    
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


