#!/bin/bash
#SBATCH --job-name=basic
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=15000

echo "Note the result may be the RC"
python ../../BASIC.py -SE ERR1421621_1M_reads.fastq.gz -b /local10G/dcroote/resources/anaconda/bin -g human -i BCR -o test
