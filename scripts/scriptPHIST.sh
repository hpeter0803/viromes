#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -J Phist
#SBATCH -o Phist.%J.out
#SBATCH -e Phist.%J.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task 28
#SBATCH --mail-user=gregoire.michoud@epfl.ch
#SBATCH --mail-type=all



source /home/${USER}/.bashrc
source activate phist

for i in finalVirus/*fna
	do mkdir ${i%_checkv_final.fna}_phist
	python /work/sber/Software/PHIST/phist.py -t28  ${i%_checkv_final.fna} MAGs/ ${i%_checkv_final.fna}_phist/common_kmers.csv ${i%_checkv_final.fna}_phist/predictions.csv
done

mkdir phist
mv finalVirus/*phist phist
