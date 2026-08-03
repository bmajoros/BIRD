#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import ProgramName

from EssexParser import EssexParser

def processEssex(filename,globalThetas):
    parser=EssexParser(filename)
    while(True):
        variant=parser.nextElem()   # returns root of the tree
        if(variant is None): break
        ID=variant.getAttribute("id")
        pools=variant.getChildren("pool")
        thetas=[]
        for pool in pools:
            poolNum=pool[0]
            poolTheta=getPoolTheta(pool)
            thetas.append([poolNum,poolTheta])
        emit(ID,thetas,globalThetas)
    parser.close()

def getPoolTheta(pool):
    dnaNodes=pool.getChildren("DNA")
    rnaNodes=pool.getChildren("RNA")
    (dnaRef,dnaAlt)=getCountSums(dnaNodes)
    (rnaRef,rnaAlt)=getCountSums(rnaNodes)
    theta=(rnaAlt/rnaRef)/(dnaAlt/dnaRef)
    return theta

def getCountSums(nodes):
    refs=[]; alts=[]
    for node in nodes:
        refs.append(int(node.getAttribute("ref")))
        alts.append(int(node.getAttribute("alt")))
    nRef=len(refs); nAlt=len(alts)
    refSum=sum(refs); altSum=sum(alts)
    normalizedRef=float(refSum+1)/float(nRef)
    normalizedAlt=float(altSum+1)/float(nAlt)
    return (normalizedRef,normalizedAlt)
        
def emit(varID,thetas,globalThetas):
    globalTheta=globalThetas[varID]
    print(varID,globalTheta,sep="\t",end="")
    for theta in thetas:
        print
    

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <input.essex> <estimates.txt>\n")
(essexFilename,estimatesFilename)=sys.argv[1:]

globalThetas=loadGlobalThetas(estimatesFilename)
processEssex(essexFilename,globalThetas)


