# snakemake rules for pulling the immcantation image depending on
# whether docker or singularity is being used


rule get_immcantation_image:
    """ pull the immcantation image using singularity

    """
    output:
        "{}/resources/repertoire/{}.sif".format(workflow.basedir, "immcantation"),
    threads: 1
    params:
        name="container_pull",
        docker_address="docker://immcantation/suite:4.1.0",
    resources:
        mem_mb=10000,
    shell:
        "module load system && "
        "mkdir -p $(dirname {output}) && cd $(dirname {output}) && "
        "img=$(basename {output}) && "
        "singularity build $img "
        "{params.docker_address}"


rule get_imgt_db:
    """ pull the imgt db

    """
    output:
        directory("{}/db/imgtdb/".format(workflow.basedir)),
    threads: 1
    params:
        name="getimgtdb",
    resources:
        mem_mb=10000,
    shell:
        "bash {params.scripts}/fetch_imgtdb.sh -o {output}"
