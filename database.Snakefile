#!/usr/bin/env python3

###########
# GLOBALS #
###########

configfile: "config.yaml"

bioconductor = ('shub://TomHarrop/singularity-containers:bioconductor_3.9'
                '@752a9788043f6a471665da4e270cd870')
biopython = ('shub://TomHarrop/singularity-containers:biopython_1.73'
             '@4a2a83e0cdff509c33227ef55906c72c')

########
# MAIN #
########

email = config['email'] if 'email' in config else ''

#########
# RULES #
#########

wildcard_constraints:
    gene_name = '\w+',
    txid = '\d+'

rule convert_fa_to_dada2:
    input:
        gb = 'output/010_database/{gene_name}-{txid}.gb',
        acc_to_taxa = 'output/010_database/{gene_name}-{txid}_acctotaxline.tab',
    output:
        fa = 'output/010_database/{gene_name}-{txid}_dada2.fa',
    log:
        'output/logs/010_database/convert_fa_to_dada2_{gene_name}-{txid}.log'
    singularity:
        biopython
    script:
        'src/convert_fa_to_dada2.py'

rule map_accession_to_taxonomy:
    input:
        nodes_dmp = 'output/010_database/ncbi-nodes-dmp.Rds',
        names_dmp = 'output/010_database/ncbi-names-dmp.Rds',
        acc_to_taxa = 'output/010_database/{gene_name}-{txid}_acctotax.csv',
    output:
        acc_to_taxline = ('output/010_database/'
                          '{gene_name}-{txid}_acctotaxline.tab')
    log:
        ('output/logs/010_database/'
         'map-accession-to-taxonomy_{gene_name}-{txid}.log')
    singularity:
        bioconductor
    script:
        'src/map_accession_to_taxonomy.R'

rule download_amp_sequences:
    params:
        search_term = lambda wildcards:
            f'{wildcards.gene_name}[Title] AND txid{wildcards.txid}[Organism]',
        email = email
    output:
        fa = 'output/010_database/{gene_name}-{txid}.fa',
        gb = 'output/010_database/{gene_name}-{txid}.gb',
        acc_to_taxa = 'output/010_database/{gene_name}-{txid}_acctotax.csv'
    log:
        'output/logs/010_database/download-amp-sequences_{gene_name}-{txid}.log'
    singularity:
        biopython
    script:
        'src/download_amp_sequences.py'

rule download_ncbi_taxonomy_db:
    output:
        nodes_dmp = 'output/010_database/ncbi-nodes-dmp.Rds',
        names_dmp = 'output/010_database/ncbi-names-dmp.Rds'
    log:
        'output/logs/010_database/download_ncbi_taxonomy_db.log'
    singularity:
        bioconductor
    script:
        'src/download_ncbi_taxonomy_db.R'
