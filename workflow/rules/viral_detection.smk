"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To identify VIRAL sequences from metagenomes
# Runs VIBRANT on metagenomes, followed by CheckV and potentially combine the both

localrules: viral_detection

rule viral_detection:
    input:
        expand(os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.{type}"), sample=SEDIMENTS, type=["faa", "fna"])
    output:
        touch("status/viral_detection.done")


#############################
# rules for viral detection #
#############################
rule vibrant:
    input:
        os.path.join(DATA_DIR, "{sample}.fa")
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
    params:
        db=config['vibrant']['db']
    wildcard_constraints:
        sample="|".join(SEDIMENTS)
    message:
        "Running VIBRANT for {wildcards.sample}"
    shell:
        "(date && export VIBRANT_DATA_PATH={params.db} && "
        "VIBRANT_run.py -t {threads} -i {input} -folder $(dirname $(dirname {output.viout1})) && date) &> {log}"

