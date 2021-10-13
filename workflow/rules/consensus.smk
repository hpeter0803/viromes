"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-13]
Latest modification:
"""

# To get consensus phage-host ids from Phist and spacepharer along with taxonomy


rule consensus:
    input:
        expand(os.path.join(RESULTS_DIR, "checkv/{sample}/{sample}_goodQual_final.fna"), sample=SAMPLES)
    output:
        touch("status/consensus.done")


##################################
# rules for Phage-host consensus #
##################################
rule phist_consensus:
    input:

    output:
        UNIQ=os.path.join(RESULTS_DIR, "consensus/{sample}/allPredictionsUniqueMatch.tsv"),
        MULT=os.path.join(RESULTS_DIR, "consensus/{sample}/allPredictionsMultipleMatch.tsv")
    conda:
        os.path.join(ENV_DIR, "renv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/consensus_phist.{sample}.log")
    message:
        ""Getting consensus Phist output for {wildcards.sample}"
    script:
        os.path.join(SRC_DIR, "getConsensusPHIST.R")

rule spacepharer_consensus:
    input:

    output:
        UNIQ=os.path.join(RESULTS_DIR, "consensus/{sample}/allPredictionsUniqueMatch.tsv"),
        MULT=os.path.join(RESULTS_DIR, "consensus/{sample}/allPredictionsMultipleMatch.tsv")
    conda:
        os.path.join(ENV_DIR, "renv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/consensus_spacepharer.{sample}.log")
    message:
        "Getting consensus spacepharer output for {wildcards.sample}"
    script:
        os.path.join(SRC_DIR, "getConsensusSpacePharer.R")

rule taxIMG_consensus:
    input:
        FILENAME, 
        IMG_HOST=os.path.join(RESULTS_DIR, "dbs/checkv-db-v1.0/IMGVR_all_Host_information.tsv"),
        IMG_SEQ=os.path.join(RESULTS_DIR, "dbs/checkv-db-v1.0/IMGVR_all_Sequence_information.tsv")
    output:

    conda:
        os.path.join(ENV_DIR, "renv.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/consensus_taxIMG.{sample}.log")
    message:
        "Getting consensus taxonomy for phage-host for {wildcards.sample}"
    script:
        os.path.join(SRC_DIR, "getConsensusTaxIMGR.R")

