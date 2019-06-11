#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=700M
#SBATCH --time=00:10:00
#SBATCH --partition=debug
#SBATCH --chdir=/scratch/stenger/SuRMISe

module load intel intel-mpi #totalview

# tvscript -event_action "error=>display_backtrace -show_arguments -show_locals" -mpi SLURM -tasks ${SLURM_NTASKS} bin/surmise -a input/clusters
# sed -i 's/^npart=.*$/npart=300000/' input/vlargeuniform.conf
# time srun -N ${SLURM_NNODES} -n ${SLURM_NTASKS} bin/surmise input/vlargeuniform
srun -N ${SLURM_NNODES} -n ${SLURM_NTASKS} amplxe-cl -collect hotspots -r /scratch/stenger/SuRMISe/vtune16/ -- bin/surmise input/clusters
