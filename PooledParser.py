#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
from EssexParser import EssexParser
from EssexNode import EssexNode
from PooledVariant import PooledVariant
from Pool import Pool
from Replicate import Replicate

#=========================================================================
# Attributes:
#   
# Instance Methods:
#   PooledParser(filename)
#   PooledVariant nextVariant() # Returns PooledVariant, or None if eof
# Class Methods:
#
# Private:
#   parser : EssexParser
#=========================================================================
class PooledParser:
    """PooledParser"""
    def __init__(self,filename):
        self.parser=EssexParser(filename)
        self.parser=EssexParser(filename)
    def nextVariant(self):
        tree=self.parser.nextElem()
        if(tree is None):
            self.parser.close()
            return None
        if(tree.getTag()!="variant"):
            raise Exception("Expecting variant")
        ID=tree.getAttribute("id")
        pooledVariant=PooledVariant(ID)
        poolSubtrees=tree.findChildren("pool")
        for subtree in poolSubtrees:
            pool=self.parsePool(subtree)
            pooledVariant.addPool(pool)
        return pooledVariant
    def parsePool(self,subtree):
        if(subtree.numElements()<4):
            raise Exception("pool contains too few fields")
        poolNumEssex=subtree[0]
        if(EssexNode.isaNode(poolNumEssex)):
            raise Exception("Expecting integer pool number as first field")
        poolNum=int(poolNumEssex)
        freq=subtree.getAttribute("freq")
        if(freq is None or freq==""):
            raise Exception("Missing 'freq' attribute for allele frequency")
        freq=float(freq)
        pool=Pool(poolNum,freq)
        dnaReps=subtree.findChildren("DNA")
        rnaReps=subtree.findChildren("RNA")
        for rep in dnaReps: pool.DNA.append(self.parseRep(rep))
        for rep in rnaReps: pool.RNA.append(self.parseRep(rep))
        return pool
    def parseRep(self,subtree):
        ref=int(subtree.getAttribute("ref")) # need to check for int
        alt=int(subtree.getAttribute("alt")) # need to check for int
        return Replicate(ref,alt)
    
