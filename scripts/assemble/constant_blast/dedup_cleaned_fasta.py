import sys
import HTSeq


def main(infile):
    """
    For some reason IMGT ships with duplicates in the fasta. This script keeps only the first entry
    """
    outfile = open(infile.split('.')[0] + '_dedup.fa', 'w')

    observed_sequences = set() 

    for entry in HTSeq.FastaReader(infile):
        if entry.name in observed_sequences:
            print(entry.name)
        else:
            entry.write_to_fasta_file(outfile)
            observed_sequences.add(entry.name)


if __name__ == '__main__':
    main(sys.argv[1])
