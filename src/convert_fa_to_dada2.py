#!/usr/bin/env python3

import csv
from Bio import SeqIO

genbank_file = snakemake.input['gb']
taxonomy_file = snakemake.input['acc_to_taxa']
fasta_file = snakemake.output['fa']

# dev
# genbank_file = 'output/gb/its1_58024.gb'
# taxonomy_file = 'output/taxonomy/taxonomy.txt'
# fasta_file = 'output/fa/its1_58024_dada2.fa'


def main():
    # read taxonomy file
    accession_to_tax = {}
    with open(taxonomy_file, 'rt') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for line in csv_reader:
            accession_to_tax[line[0]] = line[1]

    # read genbank file
    gb_records = list(SeqIO.parse(genbank_file, 'gb'))

    # update ids to match taxonomy
    for record in gb_records:
        record.description = ''
        record.id = accession_to_tax[record.id]

    # write output
    SeqIO.write(gb_records, fasta_file, 'fasta')


if __name__ == '__main__':
    main()
