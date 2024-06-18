#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
from Pool import Pool
from Replicate import Replicate
import copy

#=========================================================================
# Attributes:
#   ID : string
#   pools : array of Pool
# Instance Methods:
#   PooledVariant(ID)
#   addPool(Pool)
#   int numPools()
#   int getMaxDnaReps()
#   int getMaxRnaReps()
#   int[] getDnaReps()
#   int[] getRnaReps()
#   int[] getFreqs()
#   newVar=collapse()
#   
#   var.forceEqualFreqs()
#   bool dropHomozygousPools() # False if all pools were dropped
#   int numHomozygousPools()
#   string print()
# Class Methods:
#   
#=========================================================================
class PooledVariant:
    """PooledVariant"""
    def __init__(self,ID):
        self.ID=ID
        self.pools=[]
    def addPool(self,pool):
        self.pools.append(pool)
    def numPools(self):
        return len(self.pools)
    def getFreqs(self):
        freqs=[pool.freq for pool in self.pools]
        return freqs
    def getDnaReps(self):
        return [len(pool.DNA) for pool in self.pools]
    def getRnaReps(self):
        return [len(pool.RNA) for pool in self.pools]
    def getMaxDnaReps(self):
        return max(self.getDnaReps())
    def getMaxRnaReps(self):
        return max(self.getRnaReps())
    def getAveFreq(self): ### Assumes equal pool sizes!
        freqs=[pool.freq for pool in self.pools]
        ave=sum(freqs)/len(freqs)
        return ave
    def collapseReps(self,reps,into):
        for rep in reps:
            into.ref+=rep.ref
            into.alt+=rep.alt
    def forceEqualFreqs(self,p,q):
        for pool in self.pools:
            pool.changeFreqAndResample(p,q)
    def collapse(self):
        newVar=PooledVariant(self.ID)
        aveFreq=self.getAveFreq()
        newPool=Pool(1,aveFreq)
        dna=Replicate(0,0); rna=Replicate(0,0)
        for pool in self.pools:
            self.collapseReps(pool.DNA,dna)
            self.collapseReps(pool.RNA,rna)
        newPool.DNA.append(dna); newPool.RNA.append(rna)
        newVar.pools.append(newPool)
        return newVar
    def duplicateFirstPool(self):
        n=len(self.pools)
        first=self.pools[0]
        self.pools=[]
        for i in range(n):
            new=copy.deepcopy(first)
            self.pools.append(new)
            new.index=i
    def collapseReplicates(self):
        for pool in self.pools: pool.collapseReplicates()
    def numHomozygousPools(self):
        n=0
        for pool in self.pools:
            if(not pool.isHetPool()): n+=1
        return n
    def dropHomozygousPools(self):
        newPools=[]
        for pool in self.pools:
            #if(pool.hasHetDnaRep()): newPools.append(pool)
            if(pool.isHetPool()): newPools.append(pool)
        self.pools=newPools
        return len(newPools)>0
    def print(self):
        text="(variant (id "+self.ID+")\n"
        for pool in self.pools:
            text+="\t(pool "+str(pool.index)+"\n"
            text+="\t\t(freq "+str(pool.freq)+")\n"
            for rep in pool.DNA:
                text+="\t\t(DNA (ref "+str(rep.ref)+") (alt "+str(rep.alt)+"))\n"
            for rep in pool.RNA:
                text+="\t\t(RNA (ref "+str(rep.ref)+") (alt "+str(rep.alt)+"))\n"
            text+="\t)\n"
        text+=")"
        return text



