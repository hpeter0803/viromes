"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-13]
Latest modification:
"""

# To run SpacePHARER on viral contigs output from CheckV and MAGs
# Purpose: to identify Phage hossts


rule spacepharer_all:
    input:
        expand(os.path.join(RESULTS_DIR, "spacepharer/{sample}/{sample}_predictions.tsv"), sample=SAMPLES)
    output:
        touch("status/spacepharer.done")


#######################################
# rules for Phage-Host identification #
#######################################
rule spacepharer_dbs:
    input:
        FNA=os.path.join(RESULTS_DIR, "checkv/{sample}/{sample}_goodQual_final.fna")
    output:
        REV=directory(os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB_rev")),
        FWD=directory(os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB"))
    conda:
        os.path.join(ENV_DIR, "spacepharer.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/spacepharer_dbs.{sample}.log")
    threads:
        config["spacepharer"]["threads"]
    params:
        TMPDIR=os.path.join(RESULTS_DIR, "tmp/{sample}")
    message:
        "Creating SpacePHARER forward and reverese DBs for {wildcards.sample}"
    shell:
        "(date && "
        "spacepharer createsetdb {input.FNA} {output.REV} {params.TMPDIR} --threads {threads} --reverse-fragments 1 --compressed 1 && "
        "spacepharer createsetdb {input.FNA} {output.FWD} {params.TMPDIR} --threads {threads} --compressed 1 && "
        "date) &> {log}"

rule spacepharer:
    input:
        os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB/")
    output:
        os.path.join(RESULTS_DIR, "spacepharer/{sample}/{sample}_predictions.tsv")
    conda:
        os.path.join(ENV_DIR, "spacepharer.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/spacepharer.{sample}.log")
    threads:
        config["spacepharer"]["threads"]
    params:
        TMPDIR=os.path.join(RESULTS_DIR, "tmp/{sample}"),
        CRISPR=os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB/crisprCas/{sample}.txt")
    message:
        "Running SpacePHARER on {wildcards.sample}"
    shell:
        "(date && "
        "spacepharer easy-predict {params.CRISPR} $(dirname {input}) {output} --threads {threads} --remove-tmp-files --fmt 0 && "
        "date) &> {log}"
