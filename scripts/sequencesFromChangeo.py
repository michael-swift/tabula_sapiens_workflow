from Bio import SeqIO
import sys
import pandas as pd

bracer_file = sys.argv[1]
out_file = sys.argv[2]

bcrdf = pd.read_csv(bracer_file, sep = '\t') 
# Don't do this now because i will go through the filtering later

#combined_bcrdf = combined_bcrdf[combined_bcrdf['FUNCTIONAL'] == True]
#combined_bcrdf = combined_bcrdf[combined_bcrdf.LOCUS == 'H']

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

LRecs = [] # list of records
#Including an identifier is very important if you want to output your SeqRecord to a file. You would normally include this when creating the object:


for sequenceid in bcrdf['SEQUENCE_ID']:
    fasta_id = sequenceid
    SeqSer = bcrdf['SEQUENCE_INPUT'][bcrdf['SEQUENCE_ID'] == sequenceid]
    sample_description = bcrdf['CELL'][bcrdf['SEQUENCE_ID'] == sequenceid].values[0]
    _seq = Seq(SeqSer.values[0]) 

    simple_seq_r = SeqRecord(_seq, id = fasta_id, description = sample_description )
    LRecs.append(simple_seq_r)

from Bio import SeqIO
SeqIO.write(LRecs, out_file, "fasta")
