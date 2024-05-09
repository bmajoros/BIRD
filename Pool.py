#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
from Replicate import Replicate

#=========================================================================
# Attributes:
#   DNA : array of Replicate
#   RNA : array of Replicate
#   index : integer (which pool)
#   freq : allele frequency in pool
# Instance Methods:
#   Pool(index,freq)
#   bool isHetPool()
#   bool hasHetDnaRep() # does it have at least one het replicate?
#   pool.changeFreqAndResample(dna_freq,rna_frep)
#   pool.collapseReplicates()
#   pool.getPoolType() # 1=het, 2=homozygous ref, 3=homozygous alt
#   bool pool.isCollapsed()
# Private methods:
#   
#=========================================================================
class Pool:
    """Pool"""
    def __init__(self,index,freq):
        self.index=index
        self.DNA=[]
        self.RNA=[]
        self.freq=freq
    def isCollapsed(self):
        return len(self.DNA)==1 and len(self.RNA)==1
    def getPoolType(self):
        if(not self.isCollapsed()):
            raise Exception("Collapse reps before calling Pool.getPoolType()")
        rep=self.DNA[0]
        if(self.isHetPool()):
            return 1 # het pool
        if(rep.ref>0): return 2
        if(rep.alt>0): return 3
        raise Exception("Both alleles zero in Pool.getPoolType()")
    def addDnaRep(self,rep):
        self.DNA.append(rep)
    def addRnaRep(self,rep):
        self.RNA.append(rep)
    def isHetPool(self):
        return self.freq>0 and self.freq<1
    def hasHetDnaRep(self):
        raise Exception("Pool.hasHetDnaRep() is deprecated")
        #for x in self.DNA:
        #    if(x.isHet()): return True
        #return False
    def collapseReps(self,reps):
        ref=sum([rep.ref for rep in reps])
        alt=sum([rep.alt for rep in reps])
        rep0=reps[0]
        rep0.ref=ref; rep0.alt=alt
        del reps[1:]
    def collapseReplicates(self):
        self.collapseReps(self.DNA)
        self.collapseReps(self.RNA)
    def changeFreqAndResample(self,p,q):
        self.freq=p
        for r in self.DNA: r.resample(p)
        for r in self.RNA: r.resample(q)
                




