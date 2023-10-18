"""
Author: Gregoire Michoud
Affiliation: EPFL
Date: [2023-10-18]
Latest modification:
"""

rule taxonomy:
    input:
        os.path.join(RESULTS_DIR, "kaiju/kaiju_db_viruses.fmi"),
        os.path.join(RESULTS_DIR, "kaiju/viral_contigs_kaiju.annotated.out")
    output:
        touch("status/kaiju.done")

rule kaiju_database:
	output:
		fmi=os.path.join(RESULTS_DIR, "kaiju/kaiju_db_viruses.fmi"),
		nodes=os.path.join(RESULTS_DIR, "kaiju/nodes.dmp"),
		names=os.path.join(RESULTS_DIR, "kaiju/names.dmp")
	conda: 
        os.path.join(ENV_DIR, "kaiju.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/kaijudb.log")
	shell:
		"(date && mkdir -p $(dirname {output.fmi}) && "
		"kaiju-makedb -s viruses -t {threads} && "
		"mv viruses/kaiju_db_viruses.fmi {output.fmi} && mv nodes.dmp {output.fmi} && mv names.dmp {output.fmi} &&"
		"date) &> {log}"
	threads:
        config["kaiju"]["threads"]

rule kaiju_annotate_viruses:
	input:
		virus_fasta=rules.cat_bins.output,
		fmi=rules.kaiju_database.output.fmi,
		nodes=rules.kaiju_database.output.nodes,
		names=rules.kaiju_database.output.names
	output:
		raw=os.path.join(RESULTS_DIR, "kaiju/viral_contigs_kaiju.out"),
		annotated=os.path.join(RESULTS_DIR, "kaiju/viral_contigs_kaiju.annotated.out")
	conda: 
        os.path.join(ENV_DIR, "kaiju.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/kaiju.log")
	threads:
        config["kaiju"]["threads"]
	shell:
		"(date && kaiju -t {input.nodes} -f {input.fmi} -i {input.virus_fasta} -o {output.raw} -z {threads} && "
		"kaiju-addTaxonNames -t {input.nodes} -n {input.names} -i {output.raw} -o {output.annotated} && "
		"date) &> {log}"
