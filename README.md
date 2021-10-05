# Pipeline for annotating viral sequences based on the IMG/VR3 database
- [IMGVR3](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html)

# Setup
## Conda

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
conda env create -f requirements.yaml -n "snakemake"
```

# Running the Pipeline
Adjust the `imgvr3_config.yaml` file with the appropriate paths
- Do the following:

```bash
# Running the Pipeline
snakemake --use-conda --cores 36 --jobs 3 -s imgvr_snakefile -rp
```

## Notes:
- Database downloaded from the IMG/VR3 [website](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=IMG_VR)
- File unzipped as follows:
```bash
gunzip -c IMGVR_all_proteins.faa.gz > IMGVR_all_proteins.faa
```
