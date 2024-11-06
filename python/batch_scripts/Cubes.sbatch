#!/bin/tcsh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J ythree_job     # Job name
#SBATCH -o ythree_job.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e ythree_job.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -p wholenode          # Submit to the 'normal' or 'development' queue
#SBATCH -t 24:00:00             # Run time (hh:m) - 1.5 hours
#SBATCH -N 1                    # Total number of nodes requested (16 cores/node)
#SBATCH -n 3                   # Total number of mpi tasks requested
#SBATCH -A phy240036
set echo
set CWD = `pwd`
#set sim = `basename $CWD |cut -c -3`
set sim = `basename $CWD`
cd /home/x-ux454321/repos/p49d_turb_sims/python
ipython data_scrub_2/batch_make_cubes.py $sim

