import sys
import os

# -----------------------------
# Configuration

species=config['species']

# Cell Ranger
#### CellRanger VDJ #####


rule cellranger_vdj:
    input:
        config["raw_data"],
    output:
        "{base}/10X/{lib}/outs/metrics_summary.csv",
        "{base}/10X/{lib}/outs/filtered_contig.fasta",
        "{base}/10X/{lib}/outs/filtered_contig_annotations.csv",
    log:
        "{base}/logs/cellranger_{lib}.log",
    params:
        name="vdj_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        ref=config["cell_ranger_ref"],
    shell:
        "cd {params.base}/10X && rm -rf {wildcards.lib} && {params.cell_ranger} vdj --id={wildcards.lib} --fastqs={input} --reference={params.ref} --sample={wildcards.lib}"

rule igblast_10X:
    input:
        "{base}/10X/{lib}/outs/filtered_contig.fasta",
    output:
        temp("{base}/10X/igblast/{lib}_igblast_airr.tsv")
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        organism=config["species"],
        IGDBDIR=config["IGDBDIR"],
        seqtype=sense_lib_type
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
        "{base}/10X/igblast/{lib}_igblast_airr.tsv",
    output:
        "{base}/10X/igblast/{lib}_igblast.airr.tsv",
    run:
        df = pd.read_table(input[0], sep="\t", index_col=None)
        df["sequence_id"] = df["sequence_id"] + '_' + wildcards.lib
        df.to_csv(output[0], sep="\t", index=False)


#### bracer from Angela ####

rule get_bracer_contigs:
    input:
        get_changeos,
    output:
        fasta="{base}/SS2/{donor}/bracer_contigs.fasta",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{base}/logs/{donor}/get_fastas.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
        scripts=os.path.join(workflow.basedir, "scripts"),
    shell:
        "python {params.scripts}/get_bracer_contigs.py {input} {output}"


rule combine_bracer_contigs:
    input:
        expand("{base}/SS2/{donor}/bracer_contigs.fasta", base=base, donor=donors),
    output:
        "{base}/SS2/combined_bracer_contigs.fasta",
    conda:
        os.path.join(workflow.basedir, "envs/bcr.yaml")
    log:
        "{base}/logs/concat_fastas.log",
    params:
        name="combine_seqs",
        partition="quake,owners",
    shell:
        "cat {input} > {output}"

rule igblast_bracer_SS2:
    input:
        rules.combine_bracer_contigs.output
    output:
        "{base}/SS2/igblast/bracer_igblast.airr.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        organism="human",
        IGDBDIR=config["IGDBDIR"],
        seqtype="Ig",
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

rule igblast_tracer:
    input:
        "{base}/SS2/bracer/",
        sif="{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    output:
        "{base}/SS2/changeo/combined_db-pass.tsv",
    log:
        "{base}/logs/combined_igblast.log",
    params:
        name="extract_seqs",
        partition="quake,owners",
    shell:
        "singularity exec -B {wildcards.base}:/data {input.sif} changeo-clone -d /data/SS2/combined_bracer_contigs.fasta -n combined -o /data/SS2/changeo -p 2"

## Combined Outputs

rule combine_igblast:
    input:
        TenXs=expand(
            "{base}/10X/igblast/{lib}_igblast.airr.tsv", lib=libs, base=base
        ),
        SS2="{base}/SS2/igblast/bracer_igblast.airr.tsv",
    output:
        "{base}/vdj/combined_igblast.airr.tsv",
    log:
        "{base}/log/combineigblast.log",
    run:
        dfs = []
        infiles = input.TenXs
        infiles.append(input.SS2)
        print(infiles)
        for i in infiles:
            df = pd.read_table(i, sep = "\t")
            dfs.append(df)
        
        combined = pd.concat(dfs)
        combined.to_csv(output[0], sep="\t", index=False, header=True)
    

rule filter_vdj:
    input:
        "{base}/vdj/combined_igblast.airr.tsv",
    output:
        "{base}/vdj/combined_vdj.tsv.gz",
    conda:
        "../envs/pacbio.yaml"
    log:
        "{base}/logs/filter_vdj.log",
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/filter_igblast_output.py "
        "{input} "
        "-outdir {wildcards.base}/vdj "
        "--verbose "
        ">{log}"

rule annotate_constant_region:
    input:
        "{base}/vdj/combined_vdj.tsv.gz",
    output:
        "{base}/vdjc/combined_vdjc.tsv.gz",
    conda:
        "../envs/pacbio.yaml"
    params:
        ighc_db=config["ighc_db"]["human"],
        scripts=config["scripts"],
    log:
        "{base}/logs/annotate_constant_region.log",
    shell:
        "python {params.scripts}/blast_constant_region.py "
        "{input} "
        "-ighc_db {params.ighc_db} "
        "-outdir {wildcards.base}/vdjc "
        "> {log}"
