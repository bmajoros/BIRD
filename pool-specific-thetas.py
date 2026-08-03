#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import ProgramName

from EssexParser import EssexParser

NUM_POOLS=20

def processEssex(filename,globalThetas):
    parser=EssexParser(filename)
    while(True):
        variant=parser.nextElem()   # returns root of the tree
        if(variant is None): break
        ID=variant.getAttribute("id")
        pools=variant.findChildren("pool")
        thetas=[]
        for pool in pools:
            poolNum=int(pool[0])
            poolTheta=getPoolTheta(pool)
            poolTheta=round(poolTheta,5)
            thetas.append([poolNum,poolTheta])
        emit(ID,thetas,globalThetas)
    parser.close()

def getPoolTheta(pool):
    dnaNodes=pool.findChildren("DNA")
    rnaNodes=pool.findChildren("RNA")
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
    poolArray=["NA"]*NUM_POOLS
    for pair in thetas:
        (index,theta)=pair
        poolArray[index-1]=str(theta)
    globalTheta=globalThetas[varID]
    print(varID,globalTheta,sep="\t",end="")
    for i in range(NUM_POOLS):
        print("\t",poolArray[i],sep="",end="")
    print()

def loadGlobalThetas(filename,MIN_P_REG):
    thetas=dict()
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            varID=fields[0]; theta=float(fields[1]); Preg=float(fields[4])
            thetas[varID]=theta if Preg>=MIN_P_REG else 1.0
    return thetas

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <input.essex> <estimates.txt> <min-Preg>\n")
(essexFilename,estimatesFilename,MIN_P_REG)=sys.argv[1:]
MIN_P_REG=float(MIN_P_REG)

globalThetas=loadGlobalThetas(estimatesFilename,MIN_P_REG)
processEssex(essexFilename,globalThetas)


