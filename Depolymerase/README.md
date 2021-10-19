# Depolymerase

The file `combined.depolymerase.fasta` was made by parsing the following papers:

- Knecht et al., 2020, Front Microb,  [10.1007/s00253-015-7247-0](https://doi.org/10.3389/fmicb.2019.02949)
- Pires et al., 2016, Appl Microbiol Biotechnol, [10.1007/s00253-015-7247-0](https://doi.org/10.1007/s00253-015-7247-0)
- Latka et al, 2019, Front Microb, [10.3389/fmicb.2019.02649](https://doi.org/10.3389/fmicb.2019.02649)
- Latka et al, 2017, Appl Microbiol Biotechnol, [10.1007/s00253-017-8224-6](https://doi.org/10.1007/s00253-017-8224-6)

## PFAM of interest

[PF12219](https://pfam.xfam.org/family/PF12219#tabview=tab0)


## Steps

- Removed from the file `combined.depolymerase.fasta` uncharacterized proteins, one that was less than 50 AA long and a few of them that were mischaracterized or didn't have homologues to the file `depolymerase_clean.fasta`
- Alignment with `mafft`:

`mafft --thread 4 --auto depolymerase_clean.fasta  > depolymerase_clean_align.fasta`

- Phylogenetic tree with MEGA-X and the neighbor joining method to give the file `depolymerases_clean.tree`

- Blast against the IMGVR protein file (v2020-10-12_5.1) using `diamond` using as database the file `depolymerase_clean.fasta`

`diamond makedb -d depolymerase_clean --in depolymerase_clean.fasta`

`diamond blastp -q IMGVR_all_proteins.faa.gz -d depolymerase_clean.dmnd -f 6 -o ../depolymerase/depol_IMG.txt -e 1e-12 -k 1 --min-score 100 --id 60 -p 28 -b5 -c1`

- Got `3092` hits