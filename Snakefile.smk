from pprint import pprint
import os
import pandas as pd
import glob
from collections import defaultdict

test = False

shell.prefix("set +euo pipefail;")


configfile: "config/config.yaml"


outdir = config["out_dir"]
donors = [
    "pilot10",
    "pilot13",
    "pilot4",
    "pilot8",
    "pilot9",
    "pilot6",
    "pilot3",
    "pilot12",
]

os.makedirs(outdir, exist_ok=True)

# Extract Sample Names from file paths
FASTQ_DIR = config["raw_data"]
fq_paths = glob.glob(FASTQ_DIR + "*5prime*")

samples = []
for path in fq_paths:
    fq = path.split("/")[-1]
    sample = fq.rsplit("_", 4)[0]
    samples.append(sample)
samples = list(set(samples))


def get_changeos(wildcards):
    parentdir = wildcards["outdir"]
    donor = wildcards["donor"]
    changeo = "{}/SS2/{}/bracer/filtered_BCR_summary/changeodb.tab".format(
        parentdir, donor
    )
    return [changeo]


# Testing
if test == True:
    samples = samples[:2]


rule all:
    input:
        expand(
            "{outdir}/10X/{sample}/outs/metrics_summary.csv",
            sample=samples,
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
        summary="{outdir}/10X/{sample}/outs/metrics_summary.csv",
        fasta="{outdir}/10X/{sample}/outs/filtered_contig.fasta",
        anno="{outdir}/10X/{sample}/outs/filtered_contig_annotations.csv",
    log:
        "{outdir}/logs/cellranger_{sample}.log",
    params:
        name="vdj_cellranger",
        partition="quake",
        outdir=config["out_dir"],
        cell_ranger=config["cell_ranger"],
        ref=config["cell_ranger_ref"],
    shell:
        "cd {params.outdir}/10X && rm -rf {wildcards.sample} && {params.cell_ranger} vdj --id={wildcards.sample} --fastqs={input} --reference={params.ref} --sample={wildcards.sample}"


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

#rule combine_tracer_contigs:
    
rule changeo_igblast_bracer:
    input:
        "{outdir}/SS2/combined_bracer_contigs.fasta",
        sif="{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    output:
        "{outdir}/SS2/changeo/combined_db-pass.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{outdir}/logs/combined_igblast.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
    shell:
        "singularity exec -B {wildcards.outdir}:/data {input.sif} changeo-igblast -s /data/SS2/combined_bracer_contigs.fasta -n combined -o /data/SS2/changeo -p 2"


rule changeo_igblast_tracer:
    input:
        "{outdir}/SS2/bracer/",
        sif="{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    output:
        "{outdir}/SS2/changeo/combined_db-pass.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{outdir}/logs/combined_igblast.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
    shell:
        "singularity exec -B {wildcards.outdir}:/data {input.sif} changeo-clone -d /data//data/SS2/combined_bracer_contigs.fasta -n combined -o /data/SS2/changeo -p 2"

rule tenX_to_changeo:
    input:
        fasta="{outdir}/10X/{sample}/outs/filtered_contig.fasta",
        anno="{outdir}/10X/{sample}/outs/filtered_contig_annotations.csv",
    output:
        "{outdir}/10X/changeo/{sample}/10X_VDJs.fastq",
    log:
        "{outdir}/log_{sample}/log.log",
    shell:
        "cd {params.outdir} && rm -rf {wildcards.sample} && {params.cell_ranger} vdj --id={wildcards.sample} --fastqs={input} --reference={params.ref} --sample={wildcards.sample}"


rule combine_tenX_TCRs:
    input:
        anno="{outdir}/10X/{sample}/outs/filtered_contig_annotations.csv",
    output:
        "{outdir}/10X/changeo/{sample}/.fastq",
    log:
        "{outdir}/log_{sample}/log.log",
    shell:
        "cd {params.outdir} && rm -rf {wildcards.sample} && {params.cell_ranger} vdj --id={wildcards.sample} --fastqs={input} --reference={params.ref} --sample={wildcards.sample}"


include: "rules/get_containers.smk"
