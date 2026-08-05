#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import gzip
import ProgramName
from Rex import Rex

rex=Rex()

def processVCF(filename,MAX_VARIANTS,NUM_DONORS):
    numProcessed=0
    with gzip.open(filename,"rt") as IN:
        for line in IN:
            if(len(line<1)): continue
            if(line[0]=='#'): continue
            
            
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <in.vcf.gz> <max-variants> <#donors>\n")
(vcfFilename,MAX_VARIANTS,NUM_DONORS)=sys.argv[1:]
MAX_VARIANTS=int(MAX_VARIANTS)
NUM_DONORS=int(NUM_DONORS)

processVCF(vcfFilename,MAX_VARIANTS,NUM_DONORS)



