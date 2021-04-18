#!/usr/bin/env python3
"""
Clean IMGT germline fasta files for IgBLAST database build
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('--gapped', action='store_true', dest='gapped',
                        default=False,
                        help='Retain IMGT gaps')

    parser.add_argument('infile',
                        help='Fasta infile')

    parser.add_argument('outfile',
                        help='Cleaned fasta outfile')

    return parser.parse_args()


def main():
    # Get input and output file names
    args = parse_args()
    in_file = args.infile
    out_file = args.outfile

    # Load sequences into memory and process them
    name_set = set()
    seq_list = list()
    for rec in SeqIO.parse(in_file, 'fasta'):
        name = rec.description.split('|')[1]
        if name not in name_set and 'J-C' not in name:
            name_set.add(name)
            if args.gapped:
                seq = SeqRecord(rec.seq.upper(),
                                id=name, name=name, description=name)
            else:
                seq = SeqRecord(rec.seq.ungap('.').upper(),
                                id=name, name=name, description=name)
            seq_list.append(seq)

    # Overwrite file
    with open(out_file, 'w') as out_handle:
        writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
        writer.write_file(seq_list)

if __name__ == "__main__":
    main()
