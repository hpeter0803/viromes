"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# Taxonomic classification of VIRAL sequences identified via VIBRANT and checkV
# Runs DIAMOND against the IMG/VR3 viral taxonomy database and merges with the coverage per sample


rule taxonomy:
    input:
        expand(os.path.join(RESULTS_DIR, "diamond/{sample}.{type}"), sample=SAMPLES, type=["daa", "tsv"]),
        expand(os.path.join(RESULTS_DIR, "taxa_cov/{sample}_IMGVR_taxonomy_coverage.txt"), sample=SAMPLES)    
    output:
        touch("status/imgvr3.done")


######################################
# rules for viral TAXONOMIC analyses #
######################################
# Making DIAMOND db
rule makedb:
    input:
        IMG=os.path.join(DB_DIR, "IMGVR_all_proteins.faa")
    output:
        DB=os.path.join(RESULTS_DIR, "IMGVR_db/IMGVR_all_proteins.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/diamond_db.log")
    conda: 
        os.path.join(ENV_DIR, "diamond.yaml")
    threads:
        config["diamond"]["threads"]
    params:
        name="IMGVR"
    message: 
        "Creating the BLAST db"
    shell:
        "(date && "
        "db={output.DB} && "
        "diamond makedb --in {input.IMG} --db ${{db%.*}} --threads {threads} && "
        "date) &> {log}"

# BLAST against IMGVR3 using DIAMOND
rule blast:
    input:
        faa=rules.vibrant.output.viout2
        db=rules.makedb.output.DB
    output:
        daa=os.path.join(RESULTS_DIR, "diamond/{sample}.daa"), 
        tsv=os.path.join(RESULTS_DIR, "diamond/{sample}.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.blast.log")
    threads:
        config["diamond"]["threads"]
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    message:
        "Running Diamond BLAST against IMGVR3 for {wildcards.sample}"
    shell:
        "(date && "
        "daa={output.daa} && "
        "diamond blastp -q {input.faa} --db {input.db} --out {output.daa} -p {threads} --outfmt 100 && "
        "diamond view --daa ${{daa%.*}} --max-target-seqs 1 -p {threads} --outfmt {params.outfmt} --out {output.tsv} && "
        "date) &> {log}"


####################################
# Viral taxonomy and gene coverage #
####################################
rule taxonomy_coverage:
    input:
        ALL_SEQ=os.path.join(DB_DIR, "IMGVR_all_Sequence_information.tsv"),
        TSV=rules.blast.output.tsv,
        COV=os.path.join(COV_DIR, "{sample}_gene_coverage.txt")
    output:
        tax_cov=os.path.join(RESULTS_DIR, "taxa_cov/{sample}_IMGVR_taxonomy_coverage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.tax.cov.log")
    message:
        "Merging the IMGVR3 taxonomy and the gene coverages with BLAST output for {wildcards.sample}"
    run:
        # Importing the IMG/VR3 sequence information
        dbase=pd.read_csv(input.ALL_SEQ, sep="\t", header=0)
        
        # Selecting specific columns only
        dbase=dbase[["## UViG", "Taxonomic classification", "Host taxonomy prediction"]]

        # Importing the DIAMOND BLAST output
        tsv=pd.read_csv(input.TSV, sep="\t", header=None)
        tsv=tsv[[0, 1]]
        tsv.columns=["Contig", "## UViG"]

        # Removing extra string from UViG IDs
        tsv['## UViG']=tsv['## UViG'].str.split('|').str[0]

        # Merging the TSV file with the Taxonomy information
        merged=pd.merge(tsv, dbase, on="## UViG")

        # Removing all strings after last "_" from the contigID and writing the file
        merged['Contig']=merged['Contig'].str.rsplit('_', n=1).str.get(0)

        # Importing the gene coverage file
        cov=pd.read_csv(input.COV, sep="\t", header=0)
        cov.columns=["Contig", "Gene", "Coverage"]
        
        # Merging coverage with merged TSV+taxonomy file
        final=pd.merge(merged, cov, on="Contig")
        final.to_csv(output.tax_cov, sep="\t", header=0, index=None)
