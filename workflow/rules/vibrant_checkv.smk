"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To run VAMB on viral output from VIBRANT

localrules: concatenate, db_checkv, prep_checkv

############
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]
READS_DIR=config["reads_dir"]


rule vibrant_checkv:
    input:
        os.path.join(RESULTS_DIR, "checkv/output/goodQual_final.fna")
    output:
        touch("status/vibrant_checkv.done")


################################
# Preparing VIBRANT assemblies #
################################
rule concatenate:
    input:
        expand(os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.fna"), sample=SEDIMENTS)
    output:    
        os.path.join(RESULTS_DIR, "vamb/concatenated_viral_seqs.fna")
    threads:
        config['vamb']['threads']
    conda:
        os.path.join(ENV_DIR, "vamb.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/vamb_concat.log")
    message:
        "Concatenating Viruses for VAMB"
    shell:
        "(date && mkdir -p $(dirname {output}) && concatenate.py --nozip {output} {input} && date) &> {log}"


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
        "(date && mkdir -P $(dirname {output}) && "
        "wget -P $(dirname $(dirname {output})) https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz && "
        "cd $(dirname $(dirname {output})) && tar -zxvf checkv-db-v1.0.tar.gz && date) &> {log}"

rule prep_checkv:
    input:
        db=os.path.join(RESULTS_DIR, "dbs/checkv-db-v1.0/README.txt"),
        FNA=rules.concatenate.output
    output:
        os.path.join(RESULTS_DIR, "checkv/input/all_checkv_input.fna")
    conda:
        os.path.join(ENV_DIR, "checkv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_prep.log")
    message:
        "Editing VIBRANT output to remove spaces"
    shell:
        "(date && perl -pe 's/ (.*)_fragment_(\d+)/\_$2 $1/g' {input.FNA} > {output} && date) &> {log}"

rule checkv:
    input:
        db=os.path.join(RESULTS_DIR, "dbs/checkv-db-v1.0/README.txt"),
        FNA=rules.prep_checkv.output
    output:
        os.path.join(RESULTS_DIR, "checkv/output/quality_summary.tsv")
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
        "(date && checkv end_to_end -d $(dirname {params.DB}) {input.FNA} $(dirname {output}) -t {threads} && date) &> {log}"

rule quality_filter:
    input:
        qual_in=rules.checkv.output
    output:
        qual_out=os.path.join(RESULTS_DIR, "checkv/output/goodQual.tsv"),
        complete=os.path.join(RESULTS_DIR, "checkv/output/conplete_contigs.tsv")
    conda:
        os.path.join(ENV_DIR, "renv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_quality.log")
    message:
        "Filtering to keep the complete and non-complete calls"
    script:
        os.path.join(SRC_DIR, "checkVOutMod.R")

rule quality_final:
    input:
        TSV=rules.quality_filter.output.qual_out,
        FNA=rules.concatenate.output
    output:
        os.path.join(RESULTS_DIR, "checkv/output/goodQual_final.fna")
    conda:
        os.path.join(ENV_DIR, "samtools.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_final.log")
    message:
        "Keeping only non-complete contigs"
    shell:
        "(date && do tail -n +2 {input.TSV} | cut -f1 | samtools faidx {input.FNA} -r - > {output} && date) &> {log}"

