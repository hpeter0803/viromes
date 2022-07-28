"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-13]
Latest modification:
"""

# To run CHECKV on viral contig output from VIBRANT

localrules: db_checkv, prep_checkv, quality_filter

rule checkv_all:
    input:
        expand(os.path.join(RESULTS_DIR, "checkv/{sample}/{sample}_goodQual_final.fna"), sample=SEDIMENTS)
    output:
        touch("status/checkv.done")


#################################
# rules for viral quality check #
#################################
rule db_checkv:
    output:
        os.path.join(RESULTS_DIR, "dbs/checkv-db-v1.0/README.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_db.log")
    message:
        "Downloading the CheckV database"
    shell:
        "(date && wget -O $(dirname $(dirname {output})) https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz && "
        "cd $(dirname $(dirname {output})) && tar -zxvf checkv-db-v1.0.tar.gz && date) &> {log}"

rule prep_checkv:
    input:
        rules.vibrant.output.viout3
    output:
        os.path.join(RESULTS_DIR, "checkv/{sample}/{sample}_checkv_input.fna")
    conda:
        os.path.join(ENV_DIR, "checkv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_prep.{sample}.log")
    message:
        "Editing VIBRANT output for {wildcards.sample} to remove spaces"
    shell:
        "(date && perl -pe 's/ (.*)_fragment_(\d+)/\_$2 $1/g' {input} > {output} && date) &> {log}"

rule checkv:
    input:
        rules.prep_checkv.output
    output:
        os.path.join(RESULTS_DIR, "checkv/{sample}/quality_summary.tsv")
    conda:
        os.path.join(ENV_DIR, "checkv.yaml")
    threads:
        config['checkv']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/checkv.{sample}.log")
    params:
        DB=rules.db_checkv.output
    message:
        "Running CheckV for {wildcards.sample}"
    shell:
        "(date && checkv end_to_end -d $(dirname {params.DB}) {input} $(dirname {output}) -t {threads} && date) &> {log}"

rule quality_filter:
    input:
        qual_in=os.path.join(RESULTS_DIR, "checkv/{sample}/quality_summary.tsv")
    output:
        qual_out=os.path.join(RESULTS_DIR, "checkv/{sample}/goodQual.tsv")
    conda:
        os.path.join(ENV_DIR, "renv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_quality.{sample}.log")
    message:
        "Filtering to keep only the good quality calls for {wildcards.sample}"
    script:
        os.path.join(SRC_DIR, "checkVOutMod.R")

rule quality_final:
    input:
        TSV=rules.quality_filter.output.qual_out,
        FNA=rules.vibrant.output.viout3
    output:
        os.path.join(RESULTS_DIR, "checkv/{sample}/{sample}_goodQual_final.fna")
    conda:
        os.path.join(ENV_DIR, "samtools.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_final.{sample}.log")
    message:
        "Keeping only the good quality contigs from {wildcards.sample}"
    shell:
        "(date && do tail -n +2 {input.TSV} | cut -f1 | samtools faidx {input.FNA} -r - > {output} && date) &> {log}"
