"""
Author: Gregoire Michoud
Affiliation: EPFL
Date: [2023-10-18]
Latest modification:
"""

rule iphop:
    input:
        os.path.join(RESULTS_DIR, "iphop_db/md5checkfile.txt"),
        os.path.join(RESULTS_DIR, "iphop/Host_prediction_to_genus_m90.csv")
    output:
        touch("status/iphop.done")

rule iphop_database:
    output:
        os.path.join(RESULTS_DIR, "iphop_db/md5checkfile.txt")
    conda: 
        os.path.join(ENV_DIR, "iphop.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/iphopdb.log")
    shell:
        "(date && iphop download --db_dir $(dirname {output})"
        "date) &> {log}"

rule iphop_run:
    input:
        virus_fasta=rules.cat_bins.output,
        db=rules.iphop_database.output,
    output:
        os.path.join(RESULTS_DIR, "iphop/Host_prediction_to_genus_m90.csv")
    conda: 
        os.path.join(ENV_DIR, "iphop.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/iphop.log")
    threads:
        config["iphop"]["threads"]
    shell:
        "(date && iphop predict --fa-file {input.virus_fasta} --db_dir  $(dirname {input.db}) -o $(dirname {output}) -t {threads} && "
        "date) &> {log}"
