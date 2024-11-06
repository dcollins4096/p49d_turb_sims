#!/bin/tcsh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J batch     # Job name
#SBATCH -o batch.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e batch.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p wholenode          # Submit to the 'normal' or 'development' queue
#SBATCH -t 12:00:00             # Run time (hh:m) - 1.5 hours
#SBATCH -N 1                    # Total number of nodes requested (16 cores/node)
#SBATCH -n 1                   # Total number of mpi tasks requested
#SBATCH --mail-user=dccollins@fsu.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A phy240036  
# The next line is required if the user has more than one project

#USAGE 
#Fix the node count and time above.
#Make sure enzo is built and in the working directory.
#Make sure the parameter field on the ibrun  line.

ipython data_scrub_2/batch_make_quan_frb.py aa
#ipython data_scrub_2/batch_make_spectra.py aa
