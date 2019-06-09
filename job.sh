#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:03:00
#SBATCH --partition=debug

module load gcc/7.3.0 mvapich2 totalview

tvscript -event_action "error=>display_backtrace -show_arguments -show_locals" -mpi SLURM -tasks ${SLURM_NTASKS} bin/surmise -a input/clusters
#time srun -N ${SLURM_NNODES} -n ${SLURM_NTASKS} /scratch/stenger/SuRMISe/bin/surmise /scratch/stenger/SuRMISe/input/clusters
