#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName

VALID_COMP=("greater","less","twotailed")

def toFloat(fields):
    array=[]
    for field in fields: array.append(float(field))
    return array

def countGreater(array,threshold):
    count=0
    for x in array:
        if(x>threshold): count+=1
    return count

def countLess(array,threshold):
    count=0
    for x in array:
        if(x<threshold): count+=1
    return count

def countTwoTailed(array,threshold,threshold2):
    count=0
    for x in array:
        if(x<threshold or x>threshold2): count+=1
    return count

def analyze(fields,greaterLess,threshold,threshold2):
    array=toFloat(fields)
    n=len(array)
    count=None
    if(greaterLess=="less"): count=countLess(array,threshold)
    elif(greaterLess=="greater"): count=countGreater(array,threshold)
    else: count=countTwoTailed(array,threshold,threshold2)
    P=float(count)/float(n)
    return P
        

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <thetas.txt> <greater|less|twotailed> <threshold>\n")
(filename,greaterLess,threshold)=sys.argv[1:]
if(greaterLess not in VALID_COMP): 
    exit("comparison must be one of greater, less, or twotailed")
threshold=float(threshold)
if(threshold<0): exit("threshold cannot be negative")
threshold2=None
if(greaterLess=="twotailed"):
    if(threshold==0.0): exit("threshold cannot be 0 for twotailed")
    threshold2=float(1.0/threshold)
    if(threshold>threshold2): (threshold,threshold2)=(threshold2,threshold)

with open(filename,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        P=analyze(fields,greaterLess,threshold,threshold2)
        print(P)


