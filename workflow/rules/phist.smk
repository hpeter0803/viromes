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
        expand(os.path.join(RESULTS_DIR, "phist/{filetype}.csv"), sample=SEDIMENTS, filetype=["common_kmers", "predictions"])
    output:
        touch("status/phist.done")


#######################################
# rules for Phage-Host identification #
#######################################
rule phist:
    input:
        FNA=os.path.join(RESULTS_DIR, "vrhyme/dereplicated_bins.fna")
    output:
        KMERS=os.path.join(RESULTS_DIR, "phist/common_kmers.csv"),
        PRED=os.path.join(RESULTS_DIR, "phist/predictions.csv")
    conda:
        os.path.join(ENV_DIR, "phist.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/phist.log")
    threads:
        config["phist"]["threads"]
    params:
        MAGS=os.path.join(MAGS_DIR),
        PHIST=os.path.join(SUBMODULES, "phist/phist.py")
    message:
        "Running PHiSt on dereplicated bins"
    shell:
        "(date && mkdir -p $(dirname {output.KMERS}) && "
#        "(date && cd $(dirname {params.PHIST}) && make && "
        "python3 {params.PHIST} -t {threads} $(dirname {input.FNA}) {params.MAGS} {output.KMERS} {output.PRED} && "
        "date) &> {log}"
