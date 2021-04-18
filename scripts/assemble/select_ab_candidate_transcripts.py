import sys
import pandas as pd
from Bio import SeqIO

__author__ = 'Derek Croote'
__date__ = '2015_10_31'


def main(assembler_fasta_file, isotype_blast_tsv, outfile):
    """

    :param assembler_fasta_file: fasta file cleaned headers e.g. 1.2_1, 5.1_1, etc...
    :param isotype_blast_tsv: first column is sequence id that corresponds to the cleaned fasta header
    :param outfile: fasta written to
    :return:
    """

    ids_to_extract = blast_parse_for_fasta_seq_ids(isotype_blast_tsv)

    assembler_fasta = SeqIO.parse(assembler_fasta_file, 'fasta')

    with open(outfile, 'w') as out:
        for seq_record in assembler_fasta:
            if seq_record.id in ids_to_extract:
                out.write('>%s\n%s\n' % (seq_record.id, seq_record.seq))


def blast_parse_for_fasta_seq_ids(isotype_blast_tsv):
    try:
        blast = pd.read_csv(isotype_blast_tsv, sep='\t', header=None)
        return blast[0].unique()
    except ValueError:
        # empty file- no isotype hits
        return []

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
