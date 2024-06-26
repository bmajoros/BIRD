#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import numpy as np

#=========================================================================
# Attributes:
#   ref : int
#   alt : int
# Instance Methods:
#   Replicate(ref,alt)
#   bool rep.isHet()
#   rep.resample(freq)
# Class Methods:
#   
#=========================================================================
class Replicate:
    """Replicate"""
    def __init__(self,ref,alt):
        self.ref=ref
        self.alt=alt
    def isHet(self):
        #return self.ref>0 and self.alt>0
        raise Exception("Replicate.isHet() is deprecated")
    def resample(self,p):
        n=self.ref+self.alt
        self.alt=np.random.binomial(n,p)
        self.ref=n-self.alt


