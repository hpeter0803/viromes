#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -J CheckV
#SBATCH -o CheckV.%J.out
#SBATCH -e CheckV.%J.err
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --cpus-per-task 28
#SBATCH --mail-user=gregoire.michoud@epfl.ch
#SBATCH --mail-type=all



## Start by *phages_combined.fna files from VIBRANT

source /home/${USER}/.bashrc
source activate checkv

#Vibrant may cut the contigs in 2 or more fragments and add fragment_ at the end.
#But if they are spaces, checkV may complains hence the small modification of the file to add a number at the end of the contig number in case of fragmentation

for i in *fna
	do perl -pe 's/ (.*)_fragment_(\d+)/\_$2 $1/g' $i > t.fa
	mv t.fa $i
done

for i in *fna
	do checkv end_to_end $i ${i%.phages_combined.fna}_checkv -d /work/sber/Databases/checkv-db-v1.0/ -t 28
done


for i in *checkv/quality*
	do Rscript checkVOutMod.R $i
done

for i in *checkv_goodQual.tsv
	do tail -n +2 $i | cut -f1 | samtools faidx ${i%_checkv_goodQual.tsv}_vibrant_phages.fna -r - > ${i%_goodQual.tsv}_final.fna
done