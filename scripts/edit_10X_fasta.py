import sys

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
###########################################################################                  

parser = argparse.ArgumentParser()
parser.add_argument('input_fasta', nargs='+',
                    help="fasta to edit")
parser.add_argument('-out_fasta', default='.')
parser.add_argument('-sample', help="samplename to add to each fasta line")

args = parser.parse_args()

fasta = args.input_fasta
out_fasta= args.out_fasta
lib = args.sample

records = []
   
for seq_record in SeqIO.parse(file, 'fasta'):
    seq_record.id = seq_record.id + '_' + lib
    records.append(seq_record)

SeqIO.write(records, out_file, "fasta")
