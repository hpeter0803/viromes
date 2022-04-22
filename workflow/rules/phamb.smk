"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To run VAMB on viral output from VIBRANT

localrules: phamb, concatenate, summarise_depth

############
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]
READS_DIR=config["reads_dir"]


rule phamb:
    input:
        expand(os.path.join(RESULTS_DIR, "mapping/concatenated_viruses_{type}.txt"), type=["depth", "paired"]),
        os.path.join(RESULTS_DIR, "vamb_output/clusters.tsv")
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
        os.path.join(RESULTS_DIR, "checkv/output/goodQual_final.fna")
    output:
        os.path.join(RESULTS_DIR, "dvf/goodQual_final.fna_gt2000bp_dvfpred.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/dvf.log")
    threads:
        config["dvf"]["threads"]
    params:
        dvf=config["dvf"]["path"]
    conda:
        os.path.join(ENV_DIR, "dvf.yaml")
    message:
        "Running DeepVirFinder"
    shell:
        "(date && python3 {params.dvf} -i contigs.fna -o DVF -l 2000 -c 1 && date) &> {log}"




# Phamb run
mkdir annotations
gunzip contigs.fna.gz
python3 /user/DeepVirFinder/dvf.py -i contigs.fna -o DVF -l 2000 -c 1
mv DVF/contigs.fna_gt2000bp_dvfpred.txt annotations/all.DVF.predictions.txt
prodigal -i contigs.fna -d genes.fna -a proteins.faa -p meta -g 11
hmmsearch --cpu {threads} -E 1.0e-05 -o output.txt --tblout annotations/all.hmmMiComplete105.tbl <micompleteDB> proteins.faa
hmmsearch --cpu {threads} -E 1.0e-05 -o output.txt --tblout annotations/all.hmmVOG.tbl <VOGDB> proteins.faa
gzip contigs.fna

# 
python mag_annotation/scripts/run_RF.py contigs.fna.gz vamb/clusters.tsv annotations resultdir

ls resultsidr
resultdir/vambbins_aggregated_annotation.txt
resultdir/vambbins_RF_predictions.txt
resultsdir/vamb_bins #Concatenated predicted viral bins - writes bins in chunks to files so there might be several! 




rule complete_bins:
    input:
        TSV=rules.quality_filter.output.complete,
        FNA=rules.concatenate.output
    output:
        directory(os.path.join(RESULTS_DIR, "complete_bins"))
    conda:
        os.path.join(ENV_DIR, "biopython.yaml")
    params:
        bins=os.path.join(RESULTS_DIR, "vamb_output/bins")
    log:
        os.path.join(RESULTS_DIR, "logs/complete_bins.log")
    message:
        "Keeping only non-complete contigs"
    run:
        data=pd.read_csv(input.TSV, sep="\t", header=0, usecolss="contig_id").values
        
        with open(snakemake.input.FNA, "r") as ifile:
            for record in SeqIO.parse(ifile, "fasta"):
                if record.id in data:
                    with open(os.path.join(snakemake.output, "%.fna" % record.id), "w") as ofile:
                        SeqIO.write(record, ofile, "fasta")


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
