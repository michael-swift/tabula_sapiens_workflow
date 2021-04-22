from pprint import pprint
import os
import pandas as pd
import glob
from collections import defaultdict

test = True

shell.prefix("set +euo pipefail;")


configfile: "config/config.yaml"


outdir = config["out_dir"]
# should be loaded as config file
donors = [
    "pilot_10",
    "pilot_13",
    "pilot_4",
    "pilot_8",
    "pilot_9",
    "pilot_6",
    "pilot_3",
    "pilot_12",
]


os.makedirs(outdir, exist_ok=True)

# Extract Sample Names from file paths
FASTQ_DIR = config["raw_data"]
fq_paths = glob.glob(FASTQ_DIR + "*5prime*")

libs = []

for path in fq_paths:
    fq = path.split("/")[-1]
    lib = fq.rsplit("_", 4)[0]
    libs.append(lib)
    libs = list(set(libs))


def get_changeos(wildcards):
    parentdir = wildcards["outdir"]
    donor = wildcards["donor"]
    changeo = "{}/SS2/{}/bracer/filtered_BCR_summary/changeodb.tab".format(
        parentdir, donor
    )
    return [changeo]


# Testing
if test == True:
    libs = libs[:2]

include: "rules/get_containers.smk"

rule all:
    input:
        "{}/10X/changeo/combined_changeo.tab".format(outdir), expand(
            "{outdir}/10X/{lib}/outs/metrics_summary.csv",
            lib=libs,
            outdir=outdir,
        ),
        "{}/SS2/changeo/combined_db-pass.tsv".format(outdir),

    params:
        name="all",
        partition="normal",
    threads: 1


#### CellRanger VDJ #####

rule cellranger_vdj:
    input:
        config["raw_data"],
    output:
        summary="{outdir}/10X/{lib}/outs/metrics_summary.csv",
        fasta="{outdir}/10X/{lib}/outs/filtered_contig.fasta",
        anno="{outdir}/10X/{lib}/outs/filtered_contig_annotations.csv",
    log:
        "{outdir}/logs/cellranger_{lib}.log",
    params:
        name="vdj_cellranger",
        partition="quake",
        outdir=config["out_dir"],
        cell_ranger=config["cell_ranger"],
        ref=config["cell_ranger_ref"],
    shell:
        "cd {params.outdir}/10X && rm -rf {wildcards.lib} && {params.cell_ranger} vdj --id={wildcards.lib} --fastqs={input} --reference={params.ref} --sample={wildcards.lib}"

rule changeo_10x:
    input:
        fasta="{outdir}/10X/{lib}/outs/filtered_contig.fasta",
        anno="{outdir}/10X/{lib}/outs/filtered_contig_annotations.csv",
        sif=rules.get_immcantation_image.output
    output:
        db="{outdir}/10X/{lib}/changeo/dboutput.tab",
        directory=directory("{outdir}/10X/{lib}/changeo/"),
    log:
        "{outdir}/log_{lib}/log.log",    
    shell:
        """
        DATA_DIR={wildcards.outdir}/10X/{wildcards.lib}
        READS=/data/outs/filtered_contig.fasta
        ANNOTATIONS=/data/outs/filtered_contig_annotations.csv
        SAMPLE_NAME={wildcards.lib}
        OUT_DIR=$DATA_DIR/changeo
        DIST=0.15
        NPROC=2 
        singularity exec -B $DATA_DIR:/data {input.sif} \
        changeo-10x -s $READS -a $ANNOTATIONS -x $DIST -n $SAMPLE_NAME \
        -o $OUT_DIR -p $NPROC"""

rule combine_changeo_10X:
    input:
        expand("{outdir}/10X/{lib}/changeo/dboutput.tab", lib=libs,
            outdir=outdir)
    output:
        "{outdir}/10X/changeo/combined_changeo.tab",
    log:
        "{outdir}/log/changeocombine.log",
    shell:
        "cat {input} {output}"

#### bracer

rule get_bracer_contigs:
    input:
        get_changeos,
    output:
        fasta="{outdir}/SS2/{donor}/bracer_contigs.fasta",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{outdir}/logs/{donor}/get_fastas.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
        scripts=os.path.join(workflow.basedir, "scripts"),
    shell:
        "python {params.scripts}/get_bracer_contigs.py {input} {output}"


rule combine_bracer_contigs:
    input:
        expand("{outdir}/SS2/{donor}/bracer_contigs.fasta", outdir=outdir, donor=donors),
    output:
        "{outdir}/SS2/combined_bracer_contigs.fasta",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{outdir}/logs/concat_fastas.log",
    params:
        name="combine_seqs",
        partition="quake,owners",
    shell:
        "cat {input} > {output}"


rule changeo_igblast_bracer:
    input:
        "{outdir}/SS2/combined_bracer_contigs.fasta",
        sif="{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    output:
        "{outdir}/SS2/changeo/combined_db-pass.tsv",
    log:
        "{outdir}/logs/combined_igblast.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
    shell:
        "singularity exec -B {wildcards.outdir}:/data {input.sif} changeo-igblast -s /data/SS2/combined_bracer_contigs.fasta -n combined -o /data/SS2/changeo -p 2"

#### tracer

rule changeo_igblast_tracer:
    input:
        "{outdir}/SS2/bracer/",
        sif="{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    output:
        "{outdir}/SS2/changeo/combined_db-pass.tsv",
    log:
        "{outdir}/logs/combined_igblast.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
    shell:
        "singularity exec -B {wildcards.outdir}:/data {input.sif} changeo-clone -d /data/SS2/combined_bracer_contigs.fasta -n combined -o /data/SS2/changeo -p 2"


