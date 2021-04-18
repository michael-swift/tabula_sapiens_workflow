#!/bin/bash

# Determining whether on sherlock or singlecell
if [[ $HOSTNAME == *"singlecell"* ]]; then
    source /local10G/dcroote/resources/anaconda/envs/py3.5/bin/activate /local10G/dcroote/resources/anaconda/envs/py3.5
    SNAKEFILE=/local10G/dcroote/singlecell/singleBcell_RNA_seq/Snakefile_assembly.py
    CONFIG=/local10G/dcroote/singlecell/singleBcell_RNA_seq/assembly/run/asperg_p1.yaml

elif [[ $HOSTNAME == *"sh"* ]]; then
    source /home/dcroote/resources/anaconda/envs/py3.4/bin/activate /home/dcroote/resources/anaconda/envs/py3.4
else
    echo "I'm lost"
    exit 1
fi


# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=logs/run_assembly.$DATETIME.log


if [ $# -eq 0 ]
  then
    # Dry run snakemake
    snakemake --snakefile $SNAKEFILE --configfile $CONFIG --keep-target-files --rerun-incomplete -n -r --quiet combined_assemblies.tsv
elif [ $1 = "unlock" ]
    then
        snakemake --snakefile $SNAKEFILE --configfile $CONFIG --unlock
else

  # Run snakemake
  nohup snakemake --snakefile $SNAKEFILE combined_assemblies.tsv --configfile $CONFIG --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={resources.mem_mb} -o {params.name}.%j.log" --keep-target-files --rerun-incomplete -j 200 -w 120 -k --resources download=2 | tee $LOGFILE &
#-restart-times 2 
fi

echo Log is
echo $LOGFILE
echo
echo "Done!!"
