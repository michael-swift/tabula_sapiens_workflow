import os
import pandas as pd
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('tmp_infile')
    parser.add_argument('outfile')

    return parser.parse_args()


def main():
    """ Combines parsed assemblies into a single table. """

    args = parse_args()

    df = pd.read_table(args.tmp_infile)

    # add samplename column based on path i.e.
    # {base}/{sample}/{assembler}/{output}.tsv
    df['samplename'] = df.path.apply(lambda x: x.split('/')[-3])

    assemblies = []
    for row in df.itertuples():

        try:
            tmpdf = pd.read_table(row.path)
        except pd.errors.EmptyDataError:
            continue

        tmpdf['SAMPLENAME'] = row.samplename

        assemblies.append(tmpdf)

    outdf = pd.concat(assemblies)

    outdf.to_csv(args.outfile, sep='\t')


if __name__ == "__main__":
    main()
