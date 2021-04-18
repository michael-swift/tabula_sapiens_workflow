from __future__ import division
import pandas as pd
import sys
from Bio.Seq import Seq


def main(infile, outfile):

    df = pd.read_csv(infile, sep='\t')

    # drop unnecessary columns
    drop_cols = ['abundance', 'boundaries', 'region_lengths', 'stop_codon', 'productive', 'aa', 'reads_per_molecule',
                 'primer']
    df.drop(drop_cols, axis=1, inplace=True)

    # generate AA cdr3
    df['CDR3_seq_AA'] = df.CDR3_seq.apply(lambda x: Seq(x).translate().__str__())

    # does the AA cdr3 have a stop codon?
    df['CDR3_stop_codon'] = df.CDR3_seq_AA.str.contains('*', regex=False)
    df.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
