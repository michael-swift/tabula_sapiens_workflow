import sys
import os

include: "get_containers.smk"

# -----------------------------
# Configuration

species=config['species']

# Cell Ranger
#### CellRanger VDJ #####


rule cellranger_vdj:
    input:
        config["raw_data"],
    output:
        "{outdir}/10X/{lib}/outs/metrics_summary.csv",
        "{outdir}/10X/{lib}/outs/filtered_contig.fasta",
        "{outdir}/10X/{lib}/outs/filtered_contig_annotations.csv",
    log:
        "{outdir}/logs/cellranger_{lib}.log",
    params:
        name="vdj_cellranger",
        outdir=config["out_dir"],
        cell_ranger=config["cell_ranger"],
        ref=config["cell_ranger_ref"],
    shell:
        "cd {params.outdir}/10X && rm -rf {wildcards.lib} && {params.cell_ranger} vdj --id={wildcards.lib} --fastqs={input} --reference={params.ref} --sample={wildcards.lib}"


rule igblast_10X:
    input:
        "{outdir}/10X/{lib}/outs/filtered_contig.fasta",
    output:
        "{outdir}/10X/igblast/{lib}_igblast_airr.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        organism="hsap",
        IGDBDIR=config["IGDBDIR"],
        seqtype="ig" if "BCR" in "{lib}" else "tr",
    shell:
        """
        ml load system libuv
        wdir=$(dirname {input[0]})
        export IGDATA={params.IGDBDIR}
        igblastn \
        -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_ig_v \
        -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_ig_d \
        -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_ig_j \
        -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
        -domain_system imgt \
        -ig_seqtype {params.seqtype} \
        -organism {params.organism} \
        -outfmt 19 \
        -query {input} \
        -out {output}
        """


rule add_lib_sequence_id:
    input:
        "{outdir}/10X/igblast/{lib}_igblast_airr.tsv",
    output:
        "{outdir}/10X/igblast/{lib}_igblast.airr.tsv",
    run:
        df = pd.DataFrame(input[0], sep="\t", index_col=None)
        df["sequence_id"] = wildcards.lib
        df.to_csv(output[0], sep="\t", index=False)


#### bracer from Angela ####


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


rule igblast_SS2:
    input:
        rules.combine_bracer_contigs.output
    output:
        "{outdir}/SS2/igblast/bracer_igblast.airr.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        organism="hsap",
        IGDBDIR=config["IGDBDIR"],
        seqtype="ig" if "BCR" in "{lib}" else "tr",
    shell:
        """
        ml load system libuv
        wdir=$(dirname {input[0]})
        export IGDATA={params.IGDBDIR}
        igblastn \
        -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_ig_v \
        -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_ig_d \
        -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_ig_j \
        -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
        -domain_system imgt \
        -ig_seqtype {params.seqtype} \
        -organism {params.organism} \
        -outfmt 19 \
        -query {input} \
        -out {output}
        """


rule combine_igblast:
    input:
        TenXs=expand(
            "{outdir}/10X/igblast/{lib}_igblast.airr.tsv", lib=libs, outdir=outdir
        ),
        SS2="{outdir}/SS2/igblast/bracer_igblast.airr.tsv",
    output:
        "{outdir}/combined_igblast.airr.tsv",
    log:
        "{outdir}/log/combineigblast.log",
    run:
        df = pd.concat(input.TenXs.append(input.SS2))
        df.to_csv(output[0], sep="\t", index=False, header=True)


#### tracer from angela ####


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



# SS2

#--------------BCR------------------------
# Done by Angela        
#------------------TCR -----------------------------
# Done by Angela 