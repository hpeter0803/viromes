"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To run VAMB on viral output from VIBRANT

localrules: db_phamb

############
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]
READS_DIR=config["reads_dir"]


rule phamb:
    input:
        os.path.join(RESULTS_DIR, "phamb_output/vamb_bins"),
        os.path.join(RESULTS_DIR, "complete_bins")
    output:
        touch("status/phamb.done")


###################
# PHAMB databases #
###################
rule db_phamb:
    output:
        vog=os.path.join(RESULTS_DIR, "dbs/phamb/AllVOG.hmm"),
        bact=os.path.join(RESULTS_DIR, "dbs/phamb/Bact105.hmm")
    log:
        os.path.join(RESULTS_DIR, "logs/checkv_db.log")
    message:
        "Downloading the VOGdb and Micomplete Bacterial HMMs for PHAMB"
    shell:
        "(date && wget -O $(dirname $(dirname {output})) http://fileshare.csb.univie.ac.at/vog/vog211/vog.hmm.tar.gz && "
        "cd $(dirname $(dirname {output})) && tar -zxvf vog.hmm.tar.gz && cat vog.hmm/*.hmm > {output.vog} && "
        "wget -O $(dirname $(dirname {output.bact})) https://bitbucket.org/evolegiolab/micomplete/src/master/micomplete/share/Bact105.hmm &&"
        "date) &> {log}"


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
        FNA=rules.prodigal.output.FNA
    output:
        out=os.path.join(RESULTS_DIR, "hmmer/output.txt"),
        vog=os.path.join(RESULTS_DIR, "annotations/all.hmmVOG.tbl"),
        bact=os.path.join(RESULTS_DIR, "annotations/all.hmmMiComplete105.tbl"),
        FNA=os.path.join(RESULTS_DIR, "annotations/goodQual_final.fna.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/hmmer.log")
    threads:
        config["hmmer"]["threads"]
    conda:
        os.path.join(ENV_DIR, "hmmer.yaml")
    params:
        vog=os.path.join(RESULTS_DIR, "dbs/phamb/AllVOG.hmm"),
        bact=os.path.join(RESULTS_DIR, "dbs/phamb/Bact105.hmm")
    message:
        "Searching for VOGs and Micompete Bact105 hmms"
    shell:
        "(date && hmmsearch --cpu {threads} -E 1.0e-05 -o {output.out} --tblout {output.bact} {params.bact} {input.FAA} && "
        "hmmsearch --cpu {threads} -E 1.0e-05 -o {output.out} --tblout {output.vog} {params.vog} {input.FAA} && "
        "gzip {input.FNA} > {output.FNA} && date) &> {log}"


############
# RF model #
############
rule phamb_RF:
    input:
        FNA=rules.hmmer.output.FNA,
        cluster=rules.vamb.output,
        vog=os.path.join(RESULTS_DIR, "annotations/all.hmmVOG.tbl")
    output:
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
        "Running phamb on non-complete contigs"
    shell:
        "(date && python {params.run_RF} {input.FNA} {input.cluster} $(dirname {input.vog}) $(dirname {output}) && date) &> {log}"


# Binning COMPLETE viral contigs
rule complete_bins:
    input:
        TSV=rules.quality_filter.output.complete,
        FNA=rules.concatenate.output
    output:
        directory(os.path.join(RESULTS_DIR, "complete_bins"))
    conda:
        os.path.join(ENV_DIR, "biopython.yaml")
    params:
        bins=os.path.join(RESULTS_DIR, "phamb_output/vamb_bins")
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
