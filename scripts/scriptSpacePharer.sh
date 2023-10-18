#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH -J SpacePharer
#SBATCH -o SpacePharer.%J.out
#SBATCH -e SpacePharer.%J.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task 14
#SBATCH --mail-user=gregoire.michoud@epfl.ch
#SBATCH --mail-type=all



source /home/${USER}/.bashrc
source activate spacepharer


cd finalVirus

for i in *_checkv_final.fna
	do srun -n1 -N1 -c 14 spacepharer createsetdb ${i%_checkv_final.fna}/* ../spacePharer/${i%_checkv_final.fna}_targetSetDB_rev ../spacePharer/tmpFolder/ --threads 14 --reverse-fragments 1 --compressed 1&
done
wait

for i in *_checkv_final.fna
	do srun -n1 -N1 -c 14 spacepharer createsetdb ${i%_checkv_final.fna}/* ../spacePharer/${i%_checkv_final.fna}_targetSetDB ../spacePharer/tmpFolder/ --threads 14 --compressed 1&
done
wait

for i in *_checkv_final.fna;
	do srun -n1 -N1 -c 14 spacepharer easy-predict ../crisprCas/*txt ../spacePharer/${i%_checkv_final.fna}_targetSetDB ../spacePharer/${i%_checkv_final.fna}_predictions.tsv ../spacePharer/tmpFolder/ --threads 14 --remove-tmp-files --fmt 0&
done
wait
