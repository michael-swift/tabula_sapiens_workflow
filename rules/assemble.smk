import sys
import os
include: "common.smk"
include: "qc.smk"
include: "get_containers.smk"
# -----------------------------
# Configuration

RESOURCES='/oak/stanford/groups/quake/mswift/resources/'
MINICONDA='/oak/stanford/groups/quake/mswift/resources/miniconda3/'
ACTIVATE_CONDA='set +u && source '+MINICONDA+'bin/activate'
BRACER_ENV=MINICONDA+'envs/bracer_scratch/'
Genome=config['STAR_Genome']

if 'SRA' in config.keys() and config['SRA']:
    include: "SRA.smk"

if receptor == 'BCR':
    igblast_receptor = 'Ig'
    receptor_two_char = 'ig'
elif receptor == 'TCR':
    igblast_receptor = 'TCR'
    receptor_two_char = 'tr'
else:
    raise ValueError('Receptor must be BCR or TCR')

species=config['species']


#--------------BCR------------------------
rule baldr:
    """ Assembles B cell heavy and light chains
        Pos args for baldr:
    """
    input:
        rules.cat_R1_v3.output, 
        rules.cat_R2_v3.output, 
        sif=rules.get_BALDR_image.output
    output:"{base}/{sample}/baldr/STAR/Log.final.out"
    params:
        name="baldr_run",
        partition=default_partition,
        species='Hsap' if config['species'] == 'human' else 'Mmus'
    resources:
        mem_mb=64000
    shell:
        "mkdir -p {wildcards.base}/{wildcards.sample}/baldr && cd {wildcards.base}/{wildcards.sample}/baldr && "
        "singularity run --bind {Genome}:/genome "
        "--bind {wildcards.base}/{wildcards.sample}:/data {input.sif} "
        "/home/tools/BALDR-master/BALDR --paired /data/R1.fastq.gz,/data/R2.fastq.gz "
        "--trimmomatic /home/tools/Trimmomatic-0.38/trimmomatic-0.38.jar --adapter /home/tools/Trimmomatic-0.38/adapters/NexteraPE-PE.fa "
        " --BALDR /home/tools/BALDR-master --trinity /home/tools/trinityrnaseq-Trinity-v2.3.2/Trinity "
        "--memory 64G --threads 32 --STAR /home/tools/STAR-2.6.0c/bin/Linux_x86_64/STAR --STAR_index /genome "
        " --igblastn  /home/tools/ncbi-igblast-1.6.1/bin/igblastn"

#rule combine_fasta:
        
#------------------TCR -----------------------------

rule pull_tracer_image:
    output:
        '{}/envs/tracer-latest.img'.format(workflow.basedir)
    params:
        name='pull_tracer',
        partition="quake,owners"
    resources:
        mem_mb=8000
    shell:
        "module load system &&"
        " cd $(dirname {output}) &&"
        " img=$(basename {output}) &&"
        " singularity pull --name $img"
        " docker://teichlab/tracer"


rule tracer:
    """ TCR assembly using singularity image
        Notes:
        * tracer assemble positional args: R1, R2, cell_name, out_dir
    """
    input:
        rules.pull_tracer_image.output,
        rules.cat_R1_v3.output,
        rules.cat_R2_v3.output
    output:
        '{base}/{sample}/tracer_{sample}/filtered_TCR_seqs/filtered_TCRs.txt'
    params:
        name="tracer",
        partition="quake,owners",
        species='Hsap' if config['species'] == 'human' else 'Mmus'
    resources:
        mem_mb=10000
    shell:
        "module load system &&"
        " wdir=$(dirname $(dirname ${input[1]})) &&" ## this is a really dumb way to mount but whatever
        " mkdir -p $wdir && cd $wdir &&"
        " singularity run -B $wdir:/data"
        " {input[0]} assemble -s {params.species} -r"
        " -c /tracer/docker_helper_files/docker_tracer.conf"
        " /data/$(basename {input[1]})"
        " /data/$(basename {input[2]})"
        " tracer_{wildcards.sample} /data && touch {output}"

rule summarise_tracer_symlink:
    """ Puts symlinks to bracer assembly folders in a single
        folder as required by 'bracer summarise'

        The Python loop is used to avoid os 'Argument list too long' error when
        using a for loop in shell:
    """
    input:
        tracers=expand('{base}/{sample}/tracer_{sample}/filtered_TCR_seqs/filtered_TCRs.txt',
                       zip,
                       base=bases,
                       sample=samples)
    output:
        temp('{processed_data}/tracer/symlink.done')
    params:
        name='summarise_symlink',
        partition=default_partition,
    resources:
        mem_mb=2000
    run:
        summarise_dir = os.path.dirname(str(output))
        if not os.path.exists(summarise_dir):
            os.makedirs(summarise_dir)
        for f in input:
            cell_path = os.path.dirname(os.path.dirname(str(f)))
            cell_name = os.path.basename(cell_path)
            shell("ln -s {cell_path} {summarise_dir}/{cell_name}")
        shell("touch {output}")

rule summarise_tracer:
    input:
        img=rules.pull_tracer_image.output, go=rules.summarise_tracer_symlink.output
        #tracers=expand('{base}/{sample}/tracer_{sample}/filtered_TCR_seqs/filtered_TCRs.txt',
        #               zip,
        #               base=samplesheet.base.values.tolist(),
        #               sample=samplesheet.samplename.values.tolist())
    output:
        '{processed_data}/tracer/filtered_TCRAB_summary/recombinants.txt'
    params:
        name='summarise',
        partition='quake,owners',
        species='Hsap' if config['species'] == 'human' else 'Mmus'
    resources:
        mem_mb=8000
    shell:
        "module load system &&"
        " summarisedir=$(dirname $(dirname {output})) &&"
        " mkdir -p $summarisedir &&"
        #" for i in {input.tracers}; do"
        #" cellfld=$(dirname $(dirname $i)) &&"
        #" ln -s $cellfld $summarisedir/$(basename $cellfld); done &&"
        "cd $summarisedir &&"
        " singularity run -B $PWD:/data"
        " {input.img} summarise -s {params.species}"
        " -c /tracer/docker_helper_files/docker_tracer.conf"
        " /data"
