#!/usr/bin/env python3

import pandas
from Bio import Seq
import multiprocessing

#############
# FUNCTIONS #
#############

# quick reverse complement with biopython
def rc(x):
    return(str(Seq.Seq(x).reverse_complement()))


def get_fastq_paths(wildcards):
    sample_number = wildcards.sample_number
    if int(sample_number) < 81:
        return ({
            'r1': (f'data/reads/CL7MN-3633-{int(sample_number):02}-0-1_'
                   f'S{sample_number}_L001_R1_001.fastq.gz'),
            'r2': (f'data/reads/CL7MN-3633-{int(sample_number):02}-0-1_'
                   f'S{sample_number}_L001_R2_001.fastq.gz')})
    if int(sample_number) > 80:
        return ({
            'r1': (f'data/reads/CT9ND-3633-{int(sample_number):03}-00-01_'
                   f'S{int(sample_number) - 80}_L001_R1_001.fastq.gz'),
            'r2': (f'data/reads/CT9ND-3633-{int(sample_number):03}-00-01_'
                   f'S{int(sample_number) - 80}_L001_R2_001.fastq.gz')})


###########
# GLOBALS #
###########

p1_sample_table = 'data/OG3633-01.xlsx'
p2_sample_table = 'data/OG3633-P2.xlsx'

its_f = 'ATGCGATACTTGGTGTGAAT'
its_r = 'GACGCTTCTCCAGACTACAAT'

# singularity
bioconductor = ('shub://TomHarrop/singularity-containers:bioconductor_3.9'
                '@752a9788043f6a471665da4e270cd870')
biopython = ('shub://TomHarrop/singularity-containers:biopython_1.73'
             '@4a2a83e0cdff509c33227ef55906c72c')
cutadapt = ('shub://TomHarrop/singularity-containers:cutadapt_2.6'
            '@4154a4f91aee91acf5d5db67820eaa7a')

########
# MAIN #
########

# read the sample table
p1_sample_data = pandas.read_excel(
    p1_sample_table,
    skiprows=list(range(0, 44)),
    nrows=80)
p2_sample_data = pandas.read_excel(
    p2_sample_table,
    skiprows=list(range(0, 44)),
    nrows=96)

p1_libs = sorted(set(p1_sample_data.to_dict()['Library ID'].values()))
p2_libs = sorted(set(p2_sample_data.to_dict()['Library ID'].values()))

all_libs = sorted(set(p1_libs + p2_libs))
# all_libs = all_libs[1:5] # to subset

all_sample_nos = [int(x.split('-')[1]) for x in all_libs]

#########
# RULES #
#########

subworkflow database:
    snakefile: 'database.Snakefile'


rule target:
    input:
        'output/030_dada2/taxa_with_counts.csv'


# run dada2 workflow
rule merge_taxa_with_counts:
    input:
        nochim = 'output/030_dada2/seqtab_nochim.Rds',
        taxa = 'output/030_dada2/taxa.Rds',
    output:
        taxa_with_counts = 'output/030_dada2/taxa_with_counts.csv'
    log:
        'output/logs/030_dada2/merge_taxa_with_counts.log'
    singularity:
        bioconductor
    script:
        'src/merge_taxa_with_counts.R'

rule run_dada2:
    input:
        r1 = expand('output/020_preprocess/s{sample_number}_cutadapt_r1.fq.gz',
                    sample_number=all_sample_nos),
        r2 = expand('output/020_preprocess/s{sample_number}_cutadapt_r2.fq.gz',
                    sample_number=all_sample_nos),
        fa = database('output/010_database/its1-58024_dada2.fa')
    output:
        seqtab = 'output/030_dada2/seqtab.Rds',
        nochim = 'output/030_dada2/seqtab_nochim.Rds',
        taxa = 'output/030_dada2/taxa.Rds',
        errF_plot = 'output/030_dada2/errF.pdf',
        errR_plot = 'output/030_dada2/errR.pdf'
    log:
        'output/logs/030_dada2/run_dada2.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bioconductor
    script:
        'src/run_dada2.R'


# use cutadapt as recommended
rule cutadapt:
    input:
        r1 = 'output/020_preprocess/s{sample_number}_filtn_r1.fq.gz',
        r2 = 'output/020_preprocess/s{sample_number}_filtn_r2.fq.gz'
    output:
        r1 = 'output/020_preprocess/s{sample_number}_cutadapt_r1.fq.gz',
        r2 = 'output/020_preprocess/s{sample_number}_cutadapt_r2.fq.gz'
    params:
        primer_f = its_f,
        primer_r = its_r,
        rcf = rc(its_f),
        rcr = rc(its_r)
    log:
        'output/logs/020_preprocess/cutadapt_{sample_number}.log'
    singularity:
        cutadapt
    shell:
        'cutadapt '
        '-g {params.primer_f} '
        '-a {params.rcr} '
        '-G {params.primer_r} '
        '-A {params.rcf} '
        '-n 2 '
        '-o {output.r1} '
        '-p {output.r2} '
        '{input.r1} '
        '{input.r2} '
        '&> {log}'


# use dada2's default trimming
rule d2_filter_and_trim:
    input:
        unpack(get_fastq_paths)
    output:
        r1 = 'output/020_preprocess/s{sample_number}_filtn_r1.fq.gz',
        r2 = 'output/020_preprocess/s{sample_number}_filtn_r2.fq.gz'
    log:
        'output/logs/020_preprocess/d2-filter-and-trim_{sample_number}.log'
    singularity:
        bioconductor
    script:
        'src/d2_filter_and_trim.R'

