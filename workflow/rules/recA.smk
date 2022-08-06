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
        "Extracting the 'recA' gene coverage from all {wildcard.sample}"
    shell:
        "(date && "
        """if ! grep '{params.gene}' {input.kegg}; then echo 'NA'; fi | awk -v myfile={wildcards.sample} -vOFS="\t" '{print myfile,$NF}' > {output.recA} && """
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


#############################
# rules for PHAMB binnning #
#############################
rule deepvirfinder:
    input:
        rules.quality_final.output
    output:
        os.path.join(RESULTS_DIR, "dvf/goodQual_final.fna_gt2000bp_dvfpred.txt"),
        os.path.join(RESULTS_DIR, "annotations/all.DVF.predictions.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/dvf.log")
    threads:
        config["dvf"]["threads"]
    params:
        dvf=config["dvf"]["path"],
        length=config["dvf"]["length"]
    conda:
        os.path.join(ENV_DIR, "dvf.yaml")
    message:
        "Running DeepVirFinder"
    shell:
        "(date && python3 {params.dvf} -i {input} -o $(dirname {output[0]}) -l {params.length} -c {threads} && "
        "cp -v {output[0]} {output[1]} && date) &> {log}"

rule prodigal:
    input:
        rules.quality_final.output
    output:
        FNA=os.path.join(RESULTS_DIR, "prodigal/goodQual_final.fna"),
        FAA=os.path.join(RESULTS_DIR, "prodigal/goodQual_final.faa")
    log:
        os.path.join(RESULTS_DIR, "logs/prodigal.log")
    threads:
        config["prodigal"]["threads"]
    conda:
        os.path.join(ENV_DIR, "prodigal.yaml")
    message:
        "Running prodigal on the contigs"
    shell:
        "(date && prodigal -i {input} -d {output.FNA} -a {output.FAA} -p meta -g 11 && date) &> {log}"

rule hmmer:
    input:
        FAA=rules.prodigal.output.FAA,
        FNA=rules.prodigal.output.FNA,
        vog=rules.db_phamb.output.vog,
        bact=rules.db_phamb.output.bact
    output:
        out=os.path.join(RESULTS_DIR, "hmmer/output.txt"),
        vog=os.path.join(RESULTS_DIR, "annotations/all.hmmVOG.tbl"),
        bact=os.path.join(RESULTS_DIR, "annotations/all.hmmMiComplete105.tbl")
    log:
        os.path.join(RESULTS_DIR, "logs/hmmer.log")
    threads:
        config["hmmer"]["threads"]
    conda:
        os.path.join(ENV_DIR, "hmmer.yaml")
    params:
        vog=os.path.join(RESULTS_DIR, "dbs/phamb/AllVOG.hmm") 
    message:
        "Searching for VOGs and Micompete Bact105 hmms"
    shell:
        "(date && hmmsearch --cpu {threads} -E 1.0e-05 -o {output.out} --tblout {output.bact} {input.bact} {input.FAA} && "
        "hmmsearch --cpu {threads} -E 1.0e-05 -o {output.out} --tblout {output.vog} {input.vog} {input.FAA} && date) &> {log}"


############
# RF model #
############
rule phamb_RF:
    input:
        orig_FNA=rules.quality_final.output,
        cluster=rules.vamb.output,
        vog=rules.db_phamb.output.vog,
        dvf=os.path.join(RESULTS_DIR, "annotations/all.DVF.predictions.txt")
    output:
        FNA=os.path.join(RESULTS_DIR, "annotations/goodQual_final.fna.gz"),
        out=os.path.join(RESULTS_DIR, "phamb_output/vambbins_aggregated_annotation.txt"),
        bins=directory(os.path.join(RESULTS_DIR, "phamb_output/vamb_bins"))
    log:
        os.path.join(RESULTS_DIR, "logs/phamb.log")
    threads:
        config["phamb"]["threads"]
    conda:
        os.path.join(ENV_DIR, "phamb.yaml")
    params:
        run_RF=config["phamb"]["run_RF"]
    message:
        "Running phamb on non-complete contigs. NOTE: have to enter submodules/workflows folder to run the correct phamb model"
    shell:
        "(date && gzip -c {input.orig_FNA} > {output.FNA} && "
        "cd {params.run_RF} && "
        "python mag_annotation/scripts/run_RF.py {output.FNA} {input.cluster} $(dirname {input.dvf}) $(dirname {output.out}) && date) &> {log}"


# Binning COMPLETE viral contigs
rule complete_bins:
    input:
        TSV=rules.quality_filter.output.complete,
        FNA=rules.concatenate.output
    output:
        directory(os.path.join(RESULTS_DIR, "phamb_output/vamb_bins_complete"))
    conda:
        os.path.join(ENV_DIR, "biopython.yaml")
    params:
#        bins=os.path.join(RESULTS_DIR, "phamb_output/vamb_bins"),
        predictions=os.path.join(RESULTS_DIR, "phamb_output/vambbins_aggregated_annotation.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/complete_bins.log")
    message:
        "Separating COMPLETE bins into individual fasta files"
    script:
        os.path.join(SRC_DIR, "complete_bins.py")
        

# rule checkm:
#     input:
#         rules.vamb.output
#     output:
#         os.path.join(RESULTS_DIR, "checkm/checkm_report.txt")
#     log:
#         os.path.join(RESULTS_DIR, "logs/checkm_vamb.log")
#     threads:
#         config["checkm"]["threads"]
#     params:
#         bins=os.path.join(RESULTS_DIR, "vamb_output/bins")
#     conda:
#         os.path.join(ENV_DIR, "checkm.yaml")
#     message:
#         "Running CheckM across all samples"
#     shell:
#         "(date && checkm lineage_wf -t {threads} -x fna {params.bins} $(dirname {output}) -f {output} && "
#         "date) &> {log}"
