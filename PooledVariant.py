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
        



