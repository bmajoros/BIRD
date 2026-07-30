#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the MIT License.
# Copyright (C)2026 William H. Majoros <bmajoros@duke.edu>
#=========================================================================
import sys
import ProgramName
from Pipe import Pipe

VCF="~/lab/kd259/BIRDbath_simulations/experimental_data/chr1.vcf.gz"
OUTDIR="sim2026"
GROUPS=[1,2,4,8,16]

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <num-variants> <theta> <num-reads> <noise> <minAF> <maxAF> <pools|reps>\n")
(NUM_VAR,THETA,READS,NOISE,MIN_AF,MAX_AF,poolsOrReps)=sys.argv[1:]

for numGroups in GROUPS:
    cmd="git/sim-heterogeneity.py "+VCF+" "+NUM_VAR+" "+MIN_AF+" "+\
        MAX_AF+" "+str(numGroups)+" "+poolsOrReps+" "+THETA+" "+READS+" "+\
        NOISE+\
        " > "+OUTDIR+"/"+poolsOrReps+str(numGroups)+"-theta+"+THETA+\
        "-reads"+READS+"-noise"+NOISE+".txt"
    print(cmd)
    out=Pipe.run(cmd)
    print(out)

    



