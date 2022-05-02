"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# Taxonomic classification of VIRAL sequences identified via VIBRANT and checkV
# Runs DIAMOND against the IMG/VR3 viral taxonomy database and merges with the coverage per viral contig


rule taxonomy:
    input:
        expand(os.path.join(RESULTS_DIR, "diamond/phamb_viruses.{type}"), type=["daa", "tsv"]),
        os.path.join(RESULTS_DIR, "taxa_cov/phamb_viruses_IMGVR_taxonomy_coverage.txt")
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
        fna=rules.cat_bins.output,
        db=rules.makedb.output.DB
    output:
        daa=os.path.join(RESULTS_DIR, "diamond/phamb_viruses.daa"), 
        tsv=os.path.join(RESULTS_DIR, "diamond/phamb_viruses.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/phamb_viruses.blast.log")
    threads:
        config["diamond"]["threads"]
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    message:
        "Running Diamond BLAST against IMGVR3 for phamb viruses"
    shell:
        "(date && "
        "daa={output.daa} && "
        "diamond blastx -q {input.fna} --db {input.db} --out {output.daa} -p {threads} --outfmt 100 && "
        "diamond view --daa ${{daa%.*}} --max-target-seqs 1 -p {threads} --outfmt {params.outfmt} --out {output.tsv} && "
        "date) &> {log}"


####################################
# Viral taxonomy and gene coverage #
####################################
rule taxonomy_coverage:
    input:
        ALL_SEQ=os.path.join(DB_DIR, "IMGVR_all_Sequence_information.tsv"),
        TSV=rules.blast.output.tsv,
        COV=rules.final_summarise_depth.output.depth
    output:
        tax_cov=os.path.join(RESULTS_DIR, "taxa_cov/phamb_viruses_IMGVR_taxonomy_coverage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/phamb_viruses.tax.cov.log")
    message:
        "Merging the IMGVR3 taxonomy and the gene coverages with BLAST output for phamb viruses"
    run:
        # Importing the IMG/VR3 sequence information
        dbase=pd.read_csv(input.ALL_SEQ, sep="\t", header=0)
        
        # Selecting specific columns only
        dbase=dbase[["## UViG", "Taxonomic classification", "Host taxonomy prediction"]]

        # Importing the DIAMOND BLAST output
        tsv=pd.read_csv(input.TSV, sep="\t", header=None)
        tsv.columns=["contigName", "## UViG", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]

        # Removing extra string from UViG IDs
        tsv['## UViG']=tsv['## UViG'].str.split('|').str[0]

        # Merging the TSV file with the Taxonomy information
        merged=pd.merge(tsv, dbase, on="## UViG", how="left")

        # Importing the gene coverage file
        cov=pd.read_csv(input.COV, sep="\t", header=0)

        # remove columns with the "bam-var" suffix from the coverage file
        cov_edited=cov.loc[:, ~cov.columns.str.contains(".bam-var", case=True)]
 
        # Merging coverage with merged TSV+taxonomy file
        final=pd.merge(merged, cov_edited, on="contigName", how="left")
        final.to_csv(output.tax_cov, sep="\t", header=True, index=None)
