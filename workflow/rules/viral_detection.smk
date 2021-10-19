"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To identify VIRAL sequences from metagenomes
# Runs VIBRANT on metagenomes, followed by CheckV and potentially combine the both


rule viral_detection:
    input:
        expand(os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.simple.{type}"), sample=SAMPLES, type=["faa", "fna"]),
        expand(os.path.join(RESULTS_DIR, "checkv/{sample}/quality_summary.tsv"), sample=SAMPLES)
    output:
        touch("status/viral_detection.done")


#############################
# rules for viral detection #
#############################
rule vibrant:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Assembly/mg.assembly.merged.fa")
    output:    
        viout1=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/{sample}.prodigal.faa"),
        viout2=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.faa"),
        viout3=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.fna")
    threads:
        config['vibrant']['threads']
    conda:
        os.path.join(ENV_DIR, "vibrant.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/vibrant.{sample}.log")
    message:
        "Running VIBRANT for {wildcards.sample}"
    shell:
        "(date && python3 ./vibrant/VIBRANT/VIBRANT_run.py -t {threads} -i {input} -folder $(dirname $(dirname {output.viout1})) && date) &> {log}"

