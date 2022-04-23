"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To run VAMB on viral output from VIBRANT

localrules: cat_bins, final_checkv, prep_checkv_phamb, quality_filter_phamb

############
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]
READS_DIR=config["reads_dir"]


rule final_checkv:
    input:
        os.path.join(RESULTS_DIR, "checkv/phamb/complete_contigs.tsv")
    output:
        touch("status/final_checkv.done")


#############################################
# Preparing final run of checkv on all bins #
#############################################
rule cat_bins:
    input:
        rules.phamb_RF.output.out
    output:
        os.path.join(RESULTS_DIR, "phamb_output/all_bins.fna")
    params:
        bins=os.path.join(RESULTS_DIR, "phamb_output/vamb_bins")
    log:
        os.path.join(RESULTS_DIR, "logs/cat_bins.log")
    message:
        "Concatenating all bins for CheckV"
    shell:
        "(date && cat {params.bins}/*.fna > {output} && date) &> {log}"

rule prep_checkv_phamb:
    input:
        rules.cat_bins.output
    output:
        os.path.join(RESULTS_DIR, "checkv/phamb/all_bins.fna")
    conda:
        os.path.join(ENV_DIR, "checkv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_prep.log")
    message:
        "Editing all bins fasta to remove spaces"
    shell:
        "(date && perl -pe 's/ (.*)_fragment_(\d+)/\_$2 $1/g' {input} > {output} && date) &> {log}"

rule checkv_phamb:
    input:
        rules.prep_checkv_phamb.output
    output:
        os.path.join(RESULTS_DIR, "checkv/phamb/quality_summary.tsv")
    conda:
        os.path.join(ENV_DIR, "checkv.yaml")
    threads:
        config['checkv']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/checkv.log")
    params:
        DB=rules.db_checkv.output
    message:
        "Running CheckV"
    shell:
        "(date && checkv end_to_end -d $(dirname {params.DB}) {input} $(dirname {output}) -t {threads} && date) &> {log}"

rule quality_filter_phamb:
    input:
        qual_in=rules.checkv_phamb.output
    output:
        qual_out=os.path.join(RESULTS_DIR, "checkv/phamb/goodQual.tsv"),
        complete=os.path.join(RESULTS_DIR, "checkv/phamb/complete_contigs.tsv")
    conda:
        os.path.join(ENV_DIR, "R.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_quality.log")
    message:
        "Filtering to keep the complete and non-complete calls"
    script:
        os.path.join(SRC_DIR, "checkVOutMod.R")

