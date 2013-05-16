#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=04:00:00
#$ -l mem=1G
#$ -t 1-3600
#$ -N nmb_toy
#$ -P CosmicShear
#$ -wd /scratch/scratch/ucabtok/130516_nmb_main/
#$ -o  /scratch/scratch/ucabtok/130516_nmb_main/output/
#$ -e  /scratch/scratch/ucabtok/130516_nmb_main/output/

WDIR=~/Scratch/130516_nmb_main/
SCP_OUT=kacprzak@star.ucl.ac.uk:/import/zupcx32/kacprzak/projects/130516_nmb_main/results/
FILENAME_CONFIG=nmb_main.real.yaml
DIR_RESULTS=results
N_OBJ=1000

# -------------- do not modify beyond this point ---------------

# calculate object number
N_TASKS=$((SGE_TASK_LAST-SGE_TASK_FIRST+1))
OBJ_NUM=$((SGE_TASK_ID*N_OBJ))

# make results dicts
mkdir $WDIR/$DIR_RESULTS

# we always work in tempdir
cd $TMPDIR

# greeting
echo $SGE_TASK_ID `date` "hpc is alive"
echo WDIR is $WDIR
echo SGE_TASK_LAST $SGE_TASK_LAST
echo SGE_TASK_FIRST $SGE_TASK_FIRST
echo SGE_TASK_ID $SGE_TASK_ID
echo JOB_NAME  $JOB_NAME 
echo $OBJ_NUM
echo $N_OBJ

# load modules
echo $SGE_TASK_ID `date` "loading modules"
module load default-modules
module unload compilers
module load compilers/gnu/4.4.0
module load python/2.6.6/gnu.4.4.0
module unload compilers/gnu/4.4.0
module load compilers/intel/11.1
module load gsl/1.14/intel
module load cfitsio/3260/intel
module load fftw/3.3.1/double/intel

# create command
echo $SGE_TASK_ID `date` "creating command"
$CMD="python nmb_main run $FILENAME_CONFIG -v 0 --obj_num $OBJ_NUM --nimages $N_OBJ"
echo $SGE_TASK_ID `date` "running command"
# $CMD
echo `date` "copying"
cp $TMPDIR/results.* $WDIR/$DIR_RESULTS/

# Edit and uncomment to copy results to home.
scp $TMPDIR/results.* $SCP_OUT/

echo `date` "exiting submission script"

