import sys
import pandas as pd


def main(infile, outfile, sample_name):

    blast_df = pd.read_csv(infile, sep='\t',
                           names=['contig_id', 'blast_id', 'length', 'ident', 'mismatches', 'gaps', 'start', 'end', 'evalue'])

    if blast_df.empty:
        out = pd.DataFrame(['None', 'None', 0.0, sample_name]).T
        out.to_csv(outfile, index=False, header=None)
        return 0

    blast_expanded = blast_df.blast_id.str.split('|', expand=True)
    blast_expanded.columns = ['accession', 'gene_allele', 'species', 'func', 'region', 'span', 'num_nt'] + \
                             list('abcdefghi')

    blast_df = pd.concat([blast_df, blast_expanded], axis=1)

    out = blast_df.groupby(['region', 'gene_allele']).evalue.mean().reset_index().sort(['region', 'evalue'])
    out['sample_name'] = sample_name

    # no header for simplifying "cat" shell command
    out.to_csv(outfile, index=False, header=None)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
