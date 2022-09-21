#!/bin/tcsh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J FRB     # Job name
#SBATCH -o FRB.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e FRB.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p normal          # Submit to the 'normal' or 'development' queue
#SBATCH -N 1                    # Total number of nodes requested (16 cores/node)
#SBATCH -n 4                  # Total number of mpi tasks requested
#SBATCH -t 24:00:00             # Run time (hh:m) - 1.5 hours
#SBATCH --tasks-per-node 64
#SBATCH --mail-user=dcollins4096@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
# The next line is required if the user has more than one project
#SBATCH -A TG-AST140008      # Allocation name to charge job against


ipython make_all_frb.py 5_half >& log_5_half &
ipython make_all_frb.py 5_1 >& log_5_1 &
ipython make_all_frb.py 5_2 >& log_5_2 &
ipython make_all_frb.py 5_3 >& log_5_3 &
wait
