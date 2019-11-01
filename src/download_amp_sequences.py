#!/usr/bin/env python3

import csv
import logging
from Bio import Entrez
from Bio import SeqIO

###########
# GLOBALS #
###########


search_term = snakemake.params['search_term']
email = snakemake.params['email']
genbank_file = snakemake.output['gb']
fasta_file = snakemake.output['fa']
acc_to_taxa = snakemake.output['acc_to_taxa']

# DEV: spermatocytes: txid58024[Organism] 
# search_term = ('its1[Title] '
#                'AND txid58024[Organism]')
# genbank_file = 'output/gb/its1_58024.gb'
# fasta_file = 'output/fa/its1_58024.fa'
# acc_to_taxa = 'output/taxonomy/genbank_results.txt'


########
# MAIN #
########

def main():
    # set up log
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=snakemake.log[0],
        level=logging.INFO)

    # identify myself to Entrez
    if email == '':
        raise ValueError(('Supply a valid email '
                          '(e.g. --config email=me@email.com)')
    Entrez.email = email

    logging.info(f'search_term: {search_term}')
    logging.info(f'      email: {email}')

    # initial search to get the number of hits and webenv
    with Entrez.esearch(db='nucleotide',
                        term=search_term,
                        usehistory='y',
                        idtype='acc') as handle:
        search_results = Entrez.read(handle)
        number_of_hits = int(search_results['Count'])
        webenv = search_results['WebEnv']
        query_key = search_results['QueryKey']

    logging.info("%i hits" % number_of_hits)

    # download and parse genes 1000 at a time
    batch_size = 1000
    with open(genbank_file, 'w') as genbank_handle:
        for start in range(0, number_of_hits, batch_size):
            end = min(number_of_hits, start+batch_size)
            logging.info("Start:\t%s\nEnd:\t%s" % (start, end))
            fetch_handle = Entrez.efetch(
                    db='nucleotide',
                    idtype='acc',
                    rettype='gb',
                    retmode='text',
                    webenv=webenv,
                    query_key=query_key,
                    retstart=start,
                    retmax=batch_size)
            gb_records = SeqIO.parse(fetch_handle, 'gb')
            SeqIO.write(gb_records, genbank_handle, 'gb')

    # read genbank file
    with open(genbank_file, 'r') as genbank_handle:
        gb_records = list(SeqIO.parse(genbank_file, 'gb'))

    # rename to ID only
    for record in gb_records:
        record.description = ''

    # write output
    SeqIO.write(gb_records, fasta_file, 'fasta')

    # wtf?
    # test_blank = [x for x in gb_records if x.id == '']

    # write long
    lines = []
    for record in gb_records:
        for tax in record.annotations['taxonomy']:
            lines.append([record.annotations['organism'], record.id, tax])
    with open(acc_to_taxa, 'w') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(['species', 'accession', 'taxon'])
        csv_writer.writerows(lines)


if __name__ == '__main__':
    main()

