"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-13]
Latest modification:
"""

# To run SpacePHARER on viral contigs output from CheckV and MAGs
# Purpose: to identify Phage hossts

MAGS=glob_wildcards(os.path.join(MAGS_DIR,"{name}.fa")).name

rule spacepharer_all:
    input:
        expand(os.path.join(RESULTS_DIR, "spacepharer/{mag}/spacepharer_predictions.tsv"), mag=MAGS)
    output:
        touch("status/spacepharer.done")


#######################################
# rules for Crisp-Cas identification #
#######################################
rule minced:
    input:
        BIN=os.path.join(MAGS_DIR,"{mag}.fa")
    output:
        CRISPRCas=os.path.join(RESULTS_DIR, "spacepharer/targetSetDB/crisprCas/{mag}.txt")
    conda:
        os.path.join(ENV_DIR, "spacepharer.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/minced.{mag}.log")
    params:
        TMPDIR=os.path.join(RESULTS_DIR, "tmp/minced")
    wildcard_constraints:
        mag="|".join(MAGS)
    message:
        "Getting cripsr-cas from MAGs"
    shell:
        "(date && "
        "minced -spacers -minNR 2 {input.BIN} {output.CRISPRCas} && "
        "date) &> {log}"


#######################################
# rules for Phage-Host identification #
#######################################
rule spacepharer_dbs:
    input:
        FNA=os.path.join(RESULTS_DIR, "vrhyme/dereplicated_bins.fna")
    output:
        REV=directory(os.path.join(RESULTS_DIR, "spacepharer/targetSetDB_rev")),
        FWD=directory(os.path.join(RESULTS_DIR, "spacepharer/targetSetDB"))
    conda:
        os.path.join(ENV_DIR, "spacepharer.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/spacepharer_dbs.log")
    threads:
        config["spacepharer"]["threads"]
    params:
        TMPDIR=os.path.join(RESULTS_DIR, "tmp/spacepharer_dbs")
    message:
        "Creating SpacePHARER forward and reverese DBs for dereplicated bins"
    shell:
        "(date && "
        "spacepharer createsetdb {input.FNA} {output.REV} {params.TMPDIR} --threads {threads} --reverse-fragments 1 --compressed 1 && "
        "spacepharer createsetdb {input.FNA} {output.FWD} {params.TMPDIR} --threads {threads} --compressed 1 && "
        "date) &> {log}"

rule spacepharer:
    input:
        DB=os.path.join(RESULTS_DIR, "spacepharer/targetSetDB/"),
        CRISPR=os.path.join(RESULTS_DIR, "spacepharer/targetSetDB/crisprCas/{mag}.txt")
    output:
        os.path.join(RESULTS_DIR, "spacepharer/{mag}/spacepharer_predictions.tsv")
    conda:
        os.path.join(ENV_DIR, "spacepharer.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/spacepharer.{mag}.log")
    threads:
        config["spacepharer"]["threads"]
    params:
        TMPDIR=os.path.join(RESULTS_DIR, "tmp/spacepharer"),
        #CRISPR=os.path.join(RESULTS_DIR, "spacepharer/targetSetDB/crisprCas/dereplicated_bins.txt")
    message:
        "Running SpacePHARER on dereplicated bins and Crispr-cas"
    shell:
        "(date && "
        "spacepharer easy-predict {input.CRISPR} $(dirname {input.DB}) {output} --threads {threads} --remove-tmp-files --fmt 0 && "
        "date) &> {log}"

