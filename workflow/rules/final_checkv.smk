"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# To run VAMB on viral output from VIBRANT

localrules: cat_bins, cat_dereplicated_bins, final_checkv, prep_checkv_phamb, quality_filter_phamb

############
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]
READS_DIR=config["reads_dir"]


rule final_checkv:
    input:
        os.path.join(RESULTS_DIR, "checkv/phamb/complete_contigs.tsv"),
        os.path.join(RESULTS_DIR, "phamb_output/phamb_viruses_depth.txt")
    output:
        touch("status/final_checkv.done")


#############################################
# Preparing final run of checkv on all bins #
#############################################
rule cat_bins:
    input:
#        rules.phamb_RF.output.out
        os.path.join(RESULTS_DIR, "phamb_output/vambbins_aggregated_annotation.txt")
    output:
        os.path.join(RESULTS_DIR, "phamb_output/all_bins.fna")
    params:
        bins=os.path.join(RESULTS_DIR, "phamb_output/vamb_bins"),
        complete=os.path.join(RESULTS_DIR, "phamb_output/vamb_bins_complete")
    log:
        os.path.join(RESULTS_DIR, "logs/cat_bins.log")
    message:
        "Concatenating all bins for CheckV"
    shell:
        "(date && cat {params.bins}/*.fna {params.complete}/*.fna > {output} && date) &> {log}"

rule dereplicate:
    input:
        rules.cat_bins.output
    output:
        os.path.join(RESULTS_DIR, "vrhyme/all_bins/vRhyme_best_bins_fasta/vRhyme_bin_1.fasta")
    conda:
        os.path.join(ENV_DIR, "vrhyme.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/dereplicate_vrhyme.log")
    threads:
        config["vrhyme"]["threads"]
    message:
        "Dereplicating viral bins with vRhyme"
    shell:
        "(date && rm -vrf $(dirname $(dirname {output})) && "
        "vRhyme -i {input} -t {threads} -o $(dirname $(dirname $(dirname {output})))/all_bins --derep_only --method longest && "
        "if [-f {output}]; then echo 'Dereplication worked'; else cp -v {input} $(dirname $(dirname $(dirname {output})))/dereplicated_bins.fna; fi && "
        "date) &> {log}"    

rule cat_dereplicated_bins:
    input:
        rules.dereplicate.output
    output:
        os.path.join(RESULTS_DIR, "vrhyme/dereplicated_bins.fna")
    log:
        os.path.join(RESULTS_DIR, "logs/cat_dereplicated_bins.log")
    params:
        os.path.join(RESULTS_DIR, "vrhyme/all_bins/vRhyme_best_bins_fasta/vRhyme_bin_1.fasta")
    message:
        "Concatenating dereplicated bins for CheckV"
    shell:
        "(date && cat $(dirname {input})/*.fasta > {output} && date) &> {log}"

rule prep_checkv_phamb:
    input:
        rules.cat_dereplicated_bins.output
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
        DB=os.path.join(RESULTS_DIR, "dbs/checkv-db-v1.0/README.txt")
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


############################
# rules for final coverage #
############################
rule final_mapping:
    input:
        r1=os.path.join(READS_DIR, "{sample}_mg.r1.preprocessed.fq"),
        r2=os.path.join(READS_DIR, "{sample}_mg.r2.preprocessed.fq"),
        asm=rules.cat_dereplicated_bins.output,
        cov=os.path.join(RESULTS_DIR, "mapping/concatenated_viruses_depth.txt") # Putting this here so the other BAM files are deleted prior to running this rule. Insufficient disk space otherwise.
    output:
        os.path.join(RESULTS_DIR, "coverm_mapping/dereplicated_bins.fna.{sample}_mg.r1.preprocessed.fq.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/mapping.{sample}.log")
    threads:
        config["coverm"]["bigmem_threads"]
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")
    message:
        "Mapping short reads w/ CoverM for {wildcards.sample}"
    shell:
        "(date && "
        "coverm make -r {input.asm} -t {threads} -o $(dirname {output}) --discard-unmapped -c {input.r1} {input.r2} && "
        "date) &> {log}"

# BAM index
rule final_mapping_index:
    input:
        os.path.join(RESULTS_DIR, "coverm_mapping/dereplicated_bins.fna.{sample}_mg.r1.preprocessed.fq.bam")
    output:
        os.path.join(RESULTS_DIR, "coverm_mapping/dereplicated_bins.fna.{sample}_mg.r1.preprocessed.fq.bam.bai")
    log:
        os.path.join(RESULTS_DIR, "logs/mapping.{sample}.index.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Creating samtools index for BAM files"
    shell:
        "(date && samtools index {iinput} && date) &> {log}"

## BWA index
#rule final_mapping_index:
#    input:
#        rules.cat_dereplicated_bins.output
#    output:
#        expand(os.path.join(RESULTS_DIR, "mapping/phamb_viruses.{ext}"), ext=BWA_IDX_EXT)
#    log:
#        os.path.join(RESULTS_DIR, "logs/phamb_mapping.bwa.index.log")
#    threads:
#        config["bwa"]["threads"]
#    params:
#        idx_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
#    conda:
#        os.path.join(ENV_DIR, "bwa.yaml")
#    message:
#        "Mapping: BWA index for assembly mapping"
#    shell:
#        "(date && bwa index {input} -p {params.idx_prefix} && date) &> {log}"
#
## Short reads
#rule final_mapping:
#    input:
#        r1=os.path.join(READS_DIR, "{sample}_mg.r1.preprocessed.fq"),
#        r2=os.path.join(READS_DIR, "{sample}_mg.r2.preprocessed.fq"),
#        idx=expand(os.path.join(RESULTS_DIR, "mapping/phamb_viruses.{ext}"), ext=BWA_IDX_EXT),
#        cov=os.path.join(RESULTS_DIR, "mapping/concatenated_viruses_depth.txt")	# Putting this here so the other BAM files are deleted prior to running this rule. Insufficient disk space otherwise.
#    output:
#        os.path.join(RESULTS_DIR, "mapping/phamb_{sample}.sorted.bam")
#    log:
#        os.path.join(RESULTS_DIR, "logs/phamb_bwa.mem.{sample}.log")
#    threads:
#        config["bwa"]["map_threads"]
#    params:
#        asm=rules.cat_dereplicated_bins.output,
#        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0],
#        bam_prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
#        chunk_size=config["samtools"]["sort"]["chunk_size"]
#    conda:
#        os.path.join(ENV_DIR, "bwa.yaml")
#    message:
#        "Mapping short reads w/ BWA for {wildcards.sample}"
#    shell:
#        "(date && "
#        "bwa mem -t {threads} {params.idx_prefix} {input.r1} {input.r2} | "
#        "samtools view -@ {threads} -SbT {params.asm} | "
#        "samtools sort -@ {threads} -m {params.chunk_size} -T {params.bam_prefix} -o {output} && "
#        "date) &> {log}"

rule final_summarise_depth:
    input:
        expand(os.path.join(RESULTS_DIR, "coverm_mapping/dereplicated_bins.fna.{sample}_mg.r1.preprocessed.fq.bam"), sample=SEDIMENTS)
#        expand(os.path.join(RESULTS_DIR, "mapping/phamb_{sample}.sorted.bam"), sample=SEDIMENTS)
    output:
        depth=os.path.join(RESULTS_DIR, "phamb_output/phamb_viruses_depth.txt"),
#        paired=os.path.join(RESULTS_DIR, "phamb_output/phamb_viruses_paired.txt")	# not produced by COVERM
    log:
        os.path.join(RESULTS_DIR, "logs/phamb_summarise_depth.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")
    message:
        "Getting coverage for all the samples"
    shell:
        "(date && coverm contig -b {input} -m trimmed_mean -t {threads} -o {output.depth} && date) &> {log}"
#        "(date && jgi_summarize_bam_contig_depths --outputDepth {output.depth} --pairedContigs {output.paired} {input} && date) &> {log}"
