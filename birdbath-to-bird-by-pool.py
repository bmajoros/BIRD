#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import ProgramName
from PooledParser import PooledParser

def openOutputFiles(NUM_POOLS,filestem):
    files=[]
    for i in range(NUM_POOLS):
        filename=filestem+str(i+1)+".txt"
        OUT=open(filename,"wt")
        files.append(OUT)
    return files

def closeOutputFiles(files):
    for f in files: f.close()

def processFile(filename,NUM_POOLS,MAX_VARIANTS,files):
    poolCounts=[0]*NUM_POOLS
    parser=PooledParser(filename)
    while(True):
        var=parser.nextVariant() # Returns PooledVariant, or None if eof
        if(var is None): break
        processVariant(var,poolCounts,MAX_VARIANTS,files)
        if(reachedMaxVariants(poolCounts,MAX_VARIANTS)): break

def reachedMaxVariants(poolCounts,MAX_VARIANTS):
    for count in poolCounts:
        if(count<MAX_VARIANTS): return False
    return True

def processVariant(var,poolCounts,MAX_VARIANTS,files):
    for pool in var.pools:
        if(poolCounts[pool.index-1]>=MAX_VARIANTS): continue
        emit(var.ID,pool,files[pool.index-1])
        poolCounts[pool.index-1]+=1

def emit(ID,pool,fh):
    fields=[ID]
    appendReps(pool.DNA,fields)
    appendReps(pool.RNA,fields)
    fields=[str(x) for x in fields]
    joined="\t".join(fields)
    print(joined,file=fh)

def appendReps(reps,fields):
    fields.append(len(reps))
    for rep in reps:
        fields.append(rep.ref)
        fields.append(rep.alt)
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <in:birdbath.txt> <num-pools> <max-variants> <out:filestem>\nOutput files will be named filestem###.txt for each pool\n")
(infile,NUM_POOLS,MAX_VARIANTS,filestem)=sys.argv[1:]
NUM_POOLS=int(NUM_POOLS)
MAX_VARIANTS=int(MAX_VARIANTS)

# Open output files
files=openOutputFiles(NUM_POOLS,filestem)

# Process input file
processFile(infile,NUM_POOLS,MAX_VARIANTS,files)

# Close output files
closeOutputFiles(files)
