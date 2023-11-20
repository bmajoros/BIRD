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
        



