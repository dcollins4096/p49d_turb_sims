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
set log_frb = $CWD/three_log_frb
set log_quan = $CWD/three_log_quan
set log_spec = $CWD/three_log_spec
touch $log_frb
touch $log_quan
touch $log_spec
ipython data_scrub_2/batch_make_quan.py    $sim & #>>& $log_quan 
sleep 10
ipython data_scrub_2/batch_make_frb.py     $sim >>& $log_frb 
#ipython data_scrub_2/batch_make_spectra.py $sim >>& $log_spec
#ipython data_scrub_2/batch_make_cubes.py $sim
wait

