#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:00:20
#SBATCH --partition=debug

module load gcc/7.3.0 mvapich2 totalview

tvscript -event_action "error=>display_backtrace -show_arguments -show_locals" -mpi SLURM -tasks ${SLURM_NTASKS} /scratch/stenger/SuRMISe/bin/surmise -a /scratch/stenger/SuRMISe/input/alluniform
