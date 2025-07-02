#!/bin/bash
 
# Set the queueing options 
#SBATCH --cpus-per-task 20 --exclusive
#SBATCH -n 1
#SBATCH -t 0:10:00
#SBATCH -A comgsi
#SBATCH -J gsi_test

#-------------------------------------------
# point to the source code for executable and
# for loading the same modules used in the run that needs to be repeated
#-------------------------------------------

SRCDIR=/scratch1/BMC/wrfruc/mhu/rrfs/v0.9.9/rrfs-workflow
cd $SRCDIR
. /apps/lmod/lmod/init/sh
module use modulefiles
module load build_hera_intel
module list

set -x
  ulimit -s unlimited
  ulimit -a
  export OMP_NUM_THREADS=1
  export OMP_STACKSIZE=300M
 
FV3_EXE=${SRCDIR}/exec/lake_surgery_hrrr.exe
#FV3_EXE=${SRCDIR}/exec/fix_soil_vegetation_pct.exe
workdir=/scratch1/BMC/wrfruc/mhu/rrfs/v0.9.9/test

cd ${workdir}

rm -f lake_surgery_hrrr.exe
cp ${FV3_EXE}  lake_surgery_hrrr.exe
srun ./lake_surgery_hrrr.exe

#rm -f fix_soil_vegetation_pct.exe
#cp ${FV3_EXE}  fix_soil_vegetation_pct.exe
#srun ./fix_soil_vegetation_pct.exe
 
exit
