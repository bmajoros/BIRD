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
        self.pools=["NA"]*20
    def allNA(self):
        for x in self.pools:
            if(x!="NA"): return False
        return True
        
def loadBird(poolsDir,MIN_P_REG):
    variants=dict()
    files=[str(f) for f in Path(poolsDir).iterdir() if f.is_file()]
    for filename in files:
        #if(not rex.find("out(\d+).txt",filename)): continue
        if(not rex.find("P(\d+).out.txt",filename)): continue
        readBird(filename,variants,int(rex[1]),MIN_P_REG)
    return variants

def readBird(filename,variants,poolID,MIN_P_REG):
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=5): continue
            (ID,theta,left,right,Preg)=fields
            if(not rex.find("^chr",ID)): continue
            #if(float(Preg)<MIN_P_REG): continue
            variant=variants.get(ID,None)
            if(variant is None):
                variant=variants[ID]=Variant(ID)
            variant.pools[poolID-1]=theta
    
def processBirdBath(filename,BIRD,MIN_P_REG):
    with open(filename) as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=5): continue
            (ID,theta,left,right,Preg)=fields
            if(not rex.find("^chr",ID)): continue
            #if(float(Preg)<MIN_P_REG): continue
            variant=BIRD.get(ID,None)
            if(variant is None): continue
            if(variant.allNA()): continue
            fields=[ID,theta]; fields.extend(variant.pools)
            fields.append(Preg)
            print("\t".join(fields))

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <birdbath.txt> <bird-pools-dir> <min-Preg>\n")
(birdbathFilename,poolsDir,MIN_P_REG)=sys.argv[1:]
MIN_P_REG=float(MIN_P_REG)

BIRD=loadBird(poolsDir,MIN_P_REG)
processBirdBath(birdbathFilename,BIRD,MIN_P_REG)


