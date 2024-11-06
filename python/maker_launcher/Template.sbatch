#!/bin/tcsh
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J {{jobname}}     # Job name
#SBATCH -o {{jobname}}.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e {{jobname}}.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p wholenode          # Submit to the 'normal' or 'development' queue
#SBATCH -t 24:00:00             # Run time (hh:m) - 1.5 hours
#SBATCH -N 4                    # Total number of nodes requested (16 cores/node)
#SBATCH -n 512                   # Total number of mpi tasks requested
#SBATCH --mail-user=dccollins@fsu.edu
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A phy240036 
# The next line is required if the user has more than one project
set echo

#USAGE 
#Fix the node count and time above.
#Make sure enzo is built and in the working directory.
#Make sure the parameter field on the ibrun  line.

set MAIN = {{enzoname}}
if ( -e NOGO ) then
    echo "TOLD TO NOT GO"
    exit -1
endif
#module load hdf5/1.8
if ( -e OutputLog ) then
    module load openmpi
    set restart = `awk '{print $3}' OutputLog |tail -1`
    mpirun ./enzo.exe -d -r $restart
else
    mpirun ./enzo.exe -d $MAIN
endif
#The code should continue running until it times out, so it shouldn't get here
#unless it finished or failed.
touch NOGO
