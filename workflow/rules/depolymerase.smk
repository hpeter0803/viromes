"""
Author: Susheel Bhanu BUSI
Affiliation: Systems Ecology group LCSB UniLU
Date: [2021-10-06]
Latest modification:
"""

# Runs DIAMOND against the Depolymerases database
localrules: depolymerases 

rule depolymerases:
    input:
        expand(os.path.join(RESULTS_DIR, "diamond/phamb_depolymerases.{type}"), type=["daa", "tsv"])
    output:
        touch("status/depolymerase.done")


##########################################
# rules for viral DEPOLYMERASES analyses #
##########################################
# Making DIAMOND db
rule makedb_depolymerases:
    input:
        IMG=os.path.join(DATA_DIR, "viromes/Depolymerase/depolymerases_clean.fasta")
    output:
        DB=os.path.join(RESULTS_DIR, "diamond/depolymerase_clean.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/diamond_depolymerase.log")
    conda: 
        os.path.join(ENV_DIR, "diamond.yaml")
    threads:
        config["diamond"]["threads"]
    message: 
        "Creating the Depolymerase BLAST db"
    shell:
        "(date && "
        "db={output.DB} && "
        "diamond makedb --in {input.IMG} --db ${{db%.*}} --threads {threads} && "
        "date) &> {log}"

# BLAST against IMGVR3 using DIAMOND
rule blast_depolymerases:
    input:
        fna=rules.cat_bins.output,
        db=rules.makedb_depolymerases.output.DB
    output:
        daa=os.path.join(RESULTS_DIR, "diamond/phamb_depolymerases.daa"), 
        tsv=os.path.join(RESULTS_DIR, "diamond/phamb_depolymerases.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/phamb_depolymerases.blast.log")
    threads:
        config["diamond"]["threads"]
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    conda:
        os.path.join(ENV_DIR, "diamond.yaml")
    message:
        "Running Diamond BLAST against Depolymerases db"
    shell:
        "(date && "
        "daa={output.daa} && "
        "diamond blastx -q {input.fna} --db {input.db} --out {output.daa} -p {threads} --outfmt 100 && "
        "diamond view --daa ${{daa%.*}} --max-target-seqs 1 -p {threads} --outfmt {params.outfmt} --out {output.tsv} && "
        "date) &> {log}"

