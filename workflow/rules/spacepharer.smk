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
# rules for Crisp-Cas identification #
#######################################
rule minced:
    input:
        MAGS=glob_wildcards(MAGS_DIR,"{sample}/run1/Binning/selected_DASTool_bins/{name}.fa").name
    output:
        CRISPRCas=os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB/crisprCas/{MAGS}.txt")
    conda:
        os.path.join(ENV_DIR, "spacepharer.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/minced.{sample}.log")
    params:
        TMPDIR=os.path.join(RESULTS_DIR, "tmp/{sample}")
    message:
        "Getting cripsr-cas from {wildcards.sample} MAGs"
    shell:
        "(date && "
        "minced -spacers -minNR 2 {input.MAGs} {output.CRISPRCas} && "
        "date) &> {log}"


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
        DB=os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB/")
        CRISPR=rules.minced.output.CRISPRCas
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
        #CRISPR=os.path.join(RESULTS_DIR, "spacepharer/{sample}/targetSetDB/crisprCas/{sample}.txt")
    message:
        "Running SpacePHARER on {wildcards.sample}"
    shell:
        "(date && "
        "spacepharer easy-predict {input.CRISPR} $(dirname {input.DB}) {output} --threads {threads} --remove-tmp-files --fmt 0 && "
        "date) &> {log}"

