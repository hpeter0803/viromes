#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks=16
#SBATCH -J IMG
#SBATCH -o IMG.%J.out
#SBATCH -e IMG.%J.err
#SBATCH --time=72:00:00
###SBATCH --mem=240G
#SBATCH --cpus-per-task 7
#SBATCH --mail-user=gregoire.michoud@epfl.ch
#SBATCH --mail-type=all

source /home/${USER}/.bashrc
source activate metagenomicsVirus

diamond makedb --in IMGVR_all_proteins.faa -d IMGVR_all_proteins

## Take only viruses that were deemed medium or good quality in CheckV

for i in */*fa
	do prodigal -a $i\a -f gff -d ${i%.fa}.ffn -o ${i%.fa}.gff -i $i&
done
wait

for i in */*faa
	do srun -n1 -N1 -c7 diamond blastp -q $i -f 6 -p 7 -o ${i%faa}txt -d ../IMGVR_all_proteins.dmnd -k1 -b8 -c1 --header&
done
wait

