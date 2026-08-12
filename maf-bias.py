#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import ProgramName

from EssexParser import EssexParser

NUM_POOLS=20

def processEssex(filename):
    parser=EssexParser(filename)
    while(True):
        variant=parser.nextElem()   # returns root of the tree
        if(variant is None): break
        ID=variant.getAttribute("id")
        pools=variant.findChildren("pool")
        for pool in pools:
            poolNum=int(pool[0])
            v=pool.getAttribute("freq")
            poolTheta=processCounts(pool,v)
    parser.close()

def processCounts(pool,v):
    nodes=pool.findChildren("DNA")
    for node in nodes:
        ref=int(node.getAttribute("ref"))
        alt=int(node.getAttribute("alt"))
        #if(ref==0 or alt==0): continue # ignore dropout
        #p=(alt+1)/(ref+alt+2)
        if(ref==0 and alt==0): continue
        p=alt/(alt+ref)
        print(v,p,sep="\t")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <input.essex>\n")
(essexFilename,)=sys.argv[1:]

processEssex(essexFilename)


