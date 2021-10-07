# Rules to setup and download specific databases or tools

# Download RGI data
rule download_rgi_db:
    output:
        archive=temp(os.path.join(DB_DIR, "rgi/card-data.tar.bz2")),
        json=os.path.join(DB_DIR, "rgi/card.json")
    log:
        "logs/setup.rgi.db.log"
    params:
        db_url=config["rgi"]["db_url"]
    message:
        "Setup: download RGI data"
    shell:
        "(date && "
        "wget -O {output.archive} {params.db_url} --no-check-certificate && "
        "tar -C $(dirname {output.archive}) -xvf {output.archive} && "
        "date) &> {log}"

# Setup RGI: load required DB
# NOTE: to make sure that the same DB is used for all targets
rule setup_rgi_db:
    input:
        os.path.join(DB_DIR, "rgi/card.json")
    output:
        "status/rgi_setup.done"
    log:
        "logs/setup.rgi.setup.log"
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Setup: load RGI DB"
    shell:
        "(rgi clean --local && rgi load --card_json {input} --local && rgi database --version --local) &> {log} && touch {output}"

