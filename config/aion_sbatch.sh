#!/bin/bash -l

##############################
# SLURM
# NOTE: used for this script only, NOT for the snakemake call below

#SBATCH -J nomis_virome
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=4-00:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch

##############################
# SNAKEMAKE

# conda env name
SMK_ENV="snakemake" # USER INPUT REQUIRED
# number of cores for snakemake
SMK_CORES=60
# number of jobs for snakemake
SMK_JOBS=20
# snakemake file
SMK_SMK="workflow/Snakefile"
# config file
SMK_CONFIG="config/aion_config.yaml" # USER INPUT REQUIRED
# slurm config file
SMK_SLURM="config/aion_slurm.yaml"
# slurm cluster call
# SMK_CLUSTER="sbatch -p {cluster.partition} -q {cluster.qos} {cluster.explicit} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time} --job-name={cluster.job-name}"
SMK_CLUSTER="sbatch -A project_gfs_metag {cluster.explicit} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time} --job-name={cluster.job-name}"


##############################
# IMP

# activate the env
conda activate ${SMK_ENV}

# unlock the pipeline
snakemake -s ${SMK_SMK} --cores ${SMK_CORES} --local-cores 1 --jobs ${SMK_JOBS} \
--configfile ${SMK_CONFIG} --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline \
--cluster-config ${SMK_SLURM} --cluster "${SMK_CLUSTER}" --rerun-incomplete --rerun-triggers mtime -rp --unlock

# run the pipeline
snakemake -s ${SMK_SMK} --cores ${SMK_CORES} --local-cores 1 --jobs ${SMK_JOBS} \
--configfile ${SMK_CONFIG} --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline \
--cluster-config ${SMK_SLURM} --cluster "${SMK_CLUSTER}" --rerun-incomplete --rerun-triggers mtime -rp -k 
