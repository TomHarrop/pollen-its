#!/usr/bin/env python3

###########
# GLOBALS #
###########

bioconductor = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'

########
# MAIN #
########

#########
# RULES #
#########

subworkflow database:
    snakefile: 'database.Snakefile'


rule target:
    input:
        database('output/010_database/its1-58024_dada2.fa')