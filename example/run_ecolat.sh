#!/bin/bash
#PBS -l nodes=3:ppn=28,walltime=04:00:00
#PBS -q short

NP=$((3*28))

source /gpfs/home/ajamieson/.bashrc
cd $HOME
shopt -s expand_aliases

#Redirect output
exec >> "/gpfs/scratch/ajamieson/ecotest/ecolot_b200_t2000.out"
exec 2>&1

# Prints date and time
echo "Started: $(date)"

time mpirun -np ${NP} /gpfs/home/ajamieson/ecolattice/src/parallel/ecolattice /gpfs/home/ajamieson/ecolattice/src/parallel/parameters.dat

echo "Finished: $(date)"

