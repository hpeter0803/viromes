"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2022-08-06]
Latest modification:
"""

# To extract 'recA' gene coverage from NOMIS KEGG counts

localrules: recA_extract

############
# Params


rule recA_all:
    input:
        os.path.join(RESULTS_DIR, "recA/recA_cov_all.txt")
    output:
        touch("status/recA.done")


##################
# recA COVERAGES #
##################
rule recA_extract:
    input:
        kegg=os.path.join(KEGG_DIR, "{sample}_mg.KEGG.counts.tsv")
    output:
        recA=os.path.join(RESULTS_DIR, "recA/{sample}_recA.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/recA/{sample}_recA_cov.log")
    params:
        gene=config["recA"]["gene"]
    message:
        "Extracting the 'recA' gene coverage from all {wildcards.sample}"
    shell:
        "(date && "
        """if ! grep '{params.gene}' {input.kegg}; then echo 'NA'; fi | awk -v myfile='{wildcards.sample}' -vOFS="\\t" '{{print myfile,$NF}}' > {output.recA} && """
        "date) &> {log}"

rule cat_recA:
    input:
        expand(os.path.join(RESULTS_DIR, "recA/{sample}_recA.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "recA/recA_cov_all.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/recA/cat_all.log")
    message:
        "Concatenating all the recA coverages from each sample"
    shell:
        "(date && cat {input} > {output} && date) &> {log}"

