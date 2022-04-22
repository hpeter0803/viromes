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


rule vamb:
    input:
        expand(os.path.join(RESULTS_DIR, "mapping/concatenated_viruses_{type}.txt"), type=["depth", "paired"]),
        os.path.join(RESULTS_DIR, "vamb_output/clusters.tsv")
    output:
        touch("status/vamb.done")


###########################
# rules for VAMB binnning #
###########################
# BWA index
rule mapping_index:
    input:
        rules.quality_final.output
    output:
        expand(os.path.join(RESULTS_DIR, "mapping/concatenated_viruses.{ext}"), ext=BWA_IDX_EXT)
    log:
        os.path.join(RESULTS_DIR, "logs/mapping.bwa.index.log")
    threads:
        config["bwa"]["threads"]
    params:
        idx_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Mapping: BWA index for assembly mapping"
    shell:
        "(date && bwa index {input} -p {params.idx_prefix} && date) &> {log}"

# Short reads
rule mapping:
    input:
        r1=os.path.join(READS_DIR, "{sample}_mg.r1.preprocessed.fq"),
        r2=os.path.join(READS_DIR, "{sample}_mg.r2.preprocessed.fq"),
        idx=expand(os.path.join(RESULTS_DIR, "mapping/concatenated_viruses.{ext}"), ext=BWA_IDX_EXT),
    output:
        os.path.join(RESULTS_DIR, "mapping/{sample}.sorted.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/bwa.mem.{sample}.log")
    threads:
        config["bwa"]["map_threads"]
    params:
        asm=rules.quality_final.output,
        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0],
        bam_prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        chunk_size=config["samtools"]["sort"]["chunk_size"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Mapping short reads w/ BWA for {wildcards.sample}"
    shell:
        "(date && "
        "bwa mem -t {threads} {params.idx_prefix} {input.r1} {input.r2} | "
        "samtools view -@ {threads} -SbT {params.asm} | "
        "samtools sort -@ {threads} -m {params.chunk_size} -T {params.bam_prefix} -o {output} && "
        "date) &> {log}"

rule summarise_depth:
    input:
        expand(os.path.join(RESULTS_DIR, "mapping/{sample}.sorted.bam"), sample=SEDIMENTS)
    output:
        depth=os.path.join(RESULTS_DIR, "mapping/concatenated_viruses_depth.txt"),
        paired=os.path.join(RESULTS_DIR, "mapping/concatenated_viruses_paired.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/summarise_depth.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    message:
        "Getting coverage for all the samples"
    shell:
        "(date && "
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} --pairedContigs {output.paired} {input} && date) &> {log}"

rule vamb:
    input:
        asm=rules.quality_final.output,
        depth=rules.summarise_depth.output.depth,
    output:
        os.path.join(RESULTS_DIR, "vamb_output/clusters.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/run_vamb.log")
    threads:
        config["vamb"]["threads"]
    params:
        contigID=config["vamb"]["contigID"],
        length=config["vamb"]["length"],
        minfasta=config["vamb"]["minfasta"]
    conda:
        os.path.join(ENV_DIR, "vamb.yaml")
    message:
        "Running VAMB across all samples"
    shell:
        "(date && rm -rf $(dirname {output}) && vamb --outdir $(dirname {output}) --fasta {input.asm} --jgi {input.depth} -o {params.contigID} -m {params.length} --minfasta {params.minfasta} && date) &> {log}"
