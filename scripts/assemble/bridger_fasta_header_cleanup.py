import sys

__author__ = 'Derek'
__date__ = "2015/10/29"


def main(bridger_fasta, outfile):
    """
    transforms fasta header from Bridger outupt format: >comp75_seq0 len=2118 path=[0:2118]
    to: >75.X_1  where the first int is the comp#, second is the seq#, and the third is always 1
    (filler for abundance)
    """
    out = open(outfile, 'w')

    with open(bridger_fasta) as f:
        for line in f:
            if line[0] == '>':
                first_int = line.strip().split('p')[1].split('_')[0]
                second_int = line.strip().split('seq')[1].split(' ')[0]
                new_header = '>%s.%s_1\n' % (first_int, second_int)
                out.write(new_header)
            else:
                out.write(line)

    out.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
