#!/bin/tcsh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J FRB_6     # Job name
#SBATCH -o FRB_6.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e FRB_6.o%j # Name of stderr output file(%j expands to jobId)
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

echo "lunch time."

#ipython data_scrub_1/make_3d_spectra.py 6_half 
foreach i (log_6_half log_6_1 log_6_2)
    touch $i
    echo "=============" >> $i
    echo "=============" >> $i
    echo "=============" >> $i
end
ipython data_scrub_1/make_3d_spectra.py 6_half >>& log_6_half &
ipython data_scrub_1/make_3d_spectra.py 6_1 >>& log_6_1 &
ipython data_scrub_1/make_3d_spectra.py 6_2 >>& log_6_2 &
wait
