#!/bin/bash
#SBATCH -N 8
#SBATCH --ntasks=8
#SBATCH -J Vibrant
#SBATCH -o Vibrant.%J.out
#SBATCH -e Vibrant.%J.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task 28
#SBATCH --mem=240G
#SBATCH --mail-user=gregoire.michoud@epfl.ch
#SBATCH --mail-type=all


## Start by with contig files


source /home/${USER}/.bashrc
source activate vibrant


for i in *fa
	do srun -n1 -N1 -c28 VIBRANT_run.py -i $i -t 28 -d /work/sber/Databases/Vibrant/databases/ -m /work/sber/Databases/Vibrant/files/ -folder ${i%.fa}_temp &
done
wait
