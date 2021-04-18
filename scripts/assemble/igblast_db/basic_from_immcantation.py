from Bio import SeqIO
import subprocess

"""
Basic requires heavy and light separated and within each,
constant and V genes separated i.e.
    hv
    hc
    lv
    lc

This script combines germline sequences downloaded from
IMGT using immcantation and creates bowtie2 indices 
within the BASIC folder
"""

def main():
    genes = {}
    genes['BCR'] = {'hv': ['IGHV'],
                    'hc': ['IGHC'],
                    'lv': ['IGKV', 'IGLV'],
                    'lc': ['IGKC', 'IGLC']}

    genes['TCR'] = {'hv': ['TRDV'],
                    'hc': ['TRDC'],
                    'lv': ['TRGV'],
                    'lc': ['TRGC']}

    fld = 'germline_fastas'
    species = ['human', 'mouse']
    receptors = ['BCR', 'TCR']
    chains = ['hv', 'hc', 'lv', 'lc']

    for s in species:
        print('Species: {}'.format(s))
        for r in receptors:
            print('\tReceptor: {}'.format(r))
            for chain in chains:

                if chain.endswith('v'):
                    imm_sub = 'vdj'
                elif chain.endswith('c'):
                    imm_sub = 'constant'
                else:
                    assert False, "Should not be here"

                # combine all necessary germline seqs
                pre = '{}/{}/{}/imgt_{}'.format(fld, s, imm_sub, s) 
                fnames = ['{}_{}.fasta'.format(pre, x) for x in genes[r][chain]]
                records = []
                outfile = '../BASIC/db/{}/{}/{}'.format(s, r, chain)
                print('\t\tWriting to {}:'.format(outfile))
                for fasta in fnames:
                    print('\t\t\t{}'.format(fasta))
                    with open(fasta) as f:
                        recs = list(SeqIO.parse(f, 'fasta'))
                        if recs:
                            records += recs
                        else:
                            print('\nError: {} is empty. Terminating.'.format(fasta))
                            return 1

                # write out
                with open(outfile, 'w') as out:
                    SeqIO.write(records, out, 'fasta')

                with open('bowtie2.log', 'a') as log:
                    with open('bowtie2.err', 'a') as err:
                        subprocess.check_call(['bowtie2-build', outfile, outfile],
                                              stdout=log,
                                              stderr=err)


if __name__ == '__main__':
    main()
