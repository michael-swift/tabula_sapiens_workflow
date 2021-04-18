from Bio import SeqIO
import glob
files = glob.glob('*BCR/outs/filtered_contig.fasta')
   print(files)
   out_file = 'Test_writeBCR.fasta'
   records = []
   for file in files:
        for seq_record in SeqIO.parse(file, 'fasta'):
            seq_record.id = seq_record.id + '_' + file.split('/')[0]
            records.append(seq_record)
SeqIO.write(records, out_file, "fasta")
