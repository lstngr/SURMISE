#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00
#SBATCH --partition=debug
#SBATCH --chdir=/scratch/stenger/SuRMISe

module load gcc/7.3.0 mvapich2 totalview

tvscript -event_action "error=>display_backtrace -show_arguments -show_locals" -mpi SLURM -tasks ${SLURM_NTASKS} bin/surmise -a input/clusters
#time srun -N ${SLURM_NNODES} -n ${SLURM_NTASKS} bin/surmise input/clusters
