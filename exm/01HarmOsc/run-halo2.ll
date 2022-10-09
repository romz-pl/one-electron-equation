#!/bin/bash
#PBS -N rschr-01
#PBS -l nodes=1:ppn=4
#PBS -l walltime=01:00:00
#PBS -l mem=20GB
#PBS -m e
#PBS -m a
#PBS -r n
#PBS -A G15-9
#PBS -q halo2
##

# cd ${HOME}/rschr-5.9/exm/01HarmOsc

mpiexec ../../bin/rschr.x -a ./param.inp -b ./pot.def > rschr.out

#
# Submit the job to the queuing system:
# qsub run-halo2.ll
#
