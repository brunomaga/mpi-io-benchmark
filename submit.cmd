#!/bin/bash -l
#SBATCH --job-name="IO_Benchmark"
#SBATCH --partition=debug
#SBATCH --nodes=64
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:3:00
#SBATCH --mail-user=bruno.magalhaes@epfl.ch
#SBATCH --mail-type=ALL
#SBATCH --output=IO_Benchmark.log
#SBATCH --error=IO_Benchmark.err
#SBATCH --account=proj16

#======START=====
module load slurm

# output job info
echo "On which node your job has been scheduled :"
echo $SLURM_JOB_NODELIST
echo "Print current shell limits :"
echo
ulimit -a
echo "===== Beginning job execution ======"
echo

#bug fix
export  BGLOCKLESSMPIO_F_TYPE=0x47504653

srun -n 64 ./a.out 

rm ./output.*
