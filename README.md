# Virome workflow - ENSEMBLE project

## Setup
### Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

Getting the repository including sub-modules
```bash
git clone --recurse-submodules git@github.com:hpeter0803/viromes.git
```

Create the main `snakemake` environment

```bash
# create venv
conda env create -f envs/requirements.yaml -n "snakemake"
```

## Running the Pipeline
Adjust the `config/config.yaml` file with the appropriate paths
- Do the following:

```bash
# Running the Pipeline
snakemake --use-conda --cores 36 --jobs 3 -s workflows/Snakefile -rp	

# (OR)

./config/sbatch.sh 	# works only on a SLURM-enabled system
```

