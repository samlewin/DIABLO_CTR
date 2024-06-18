#!/bin/bash
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-ccl  
# module load rhel7/default-ccl            # REQUIRED - loads the basic environment
module load fftw2/intel/64/double/2.1.5
module load hdf5/impi/1.8.16
ulimit -s unlimited      