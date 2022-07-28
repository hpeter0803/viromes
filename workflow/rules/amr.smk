"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-07]
Latest modification:
"""

# To identify AMR in viral contigs
# Runs RGItool against the CARD database


rule AMR:
    input:
        expand(os.path.join(RESULTS_DIR, "rgi/{sample}/{sample}_rgi.txt"), sample=SEDIMENTS)
    output:
        touch("status/AMR.done")


###########################
# rules for AMR detection #
###########################
rule rgi:
    input:
        faa=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.faa"),
        db=os.path.join(DB_DIR, "rgi/card.json"),
        setup="status/rgi_setup.done"	# NOTE: to make sure that the same DB is used for all targets
    output:    
        faa=temp(os.path.join(RESULTS_DIR, "rgi/{sample}/input.faa")),
        txt=os.path.join(RESULTS_DIR, "rgi/{sample}/{sample}_rgi.txt")
    threads:
        config['rgi']['threads']
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    params:
        alignment_tool="DIAMOND", 
        obname=lambda wildcards, output: os.path.splitext(output.txt)[0]
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.{sample}.log")
    message:
        "Annotating for AMR in {wildcards.sample}"
    shell:
        "(date && "
        "sed 's/\*$//' {input.faa} > {output.faa} && "
        "rgi database --version --local && "
        # NOTE: https://github.com/arpcard/rgi/issues/93: KeyError: 'snp' --> re-run
        "rgi main --input_sequence {output.faa} --output_file {params.obname} --input_type protein --local -a {params.alignment_tool} --clean -n {threads} || "
        "rgi main --input_sequence {output.faa} --output_file {params.obname} --input_type protein --local -a {params.alignment_tool} --clean -n {threads} && "
        "date) &> {log}"
