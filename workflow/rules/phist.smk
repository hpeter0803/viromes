"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-13]
Latest modification:
"""

# To run PHIST on viral contigs output from CheckV and MAGs
# Purpose: to identify Phage hossts


rule phist_all:
    input:
        expand(os.path.join(RESULTS_DIR, "phist/{sample}/{filetype}.csv"), sample=SEDIMENTS, filetype=["common_kmers", "predictions"])
    output:
        touch("status/phist.done")


#######################################
# rules for Phage-Host identification #
#######################################
rule phist:
    input:
        FNA=os.path.join(RESULTS_DIR, "checkv/{sample}/{sample}_goodQual_final.fna")
    output:
        KMERS=os.path.join(RESULTS_DIR, "phist/{sample}/common_kmers.csv"),
        PRED=os.path.join(RESULTS_DIR, "phist/{sample}/predictions.csv")
    conda:
        os.path.join(ENV_DIR, "phist.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/phist.{sample}.log")
    threads:
        config["phist"]["threads"]
    params:
        MAGS=os.path.join(MAGS_DIR)
    message:
        "Running PHiSt on {wildcards.sample}"
    shell:
        "(date && mkdir $(dirname {output.KMERS} && "
        "python phist.py -t {threads} {input.FNA} $(dirname {params.MAGS}) {output.KMERS} {output.PRED} && "
        "date) &> {log}"
