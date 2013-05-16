#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:10:0
#$ -l mem=1G
#$ -t 1-4
#$ -N nmb_main
#$ -P CosmicShear
#$ -wd /scratch/scratch/ucabtok/130515_nmb_main/000/
#$ -o  /scratch/scratch/ucabtok/130515_nmb_main/000/output/
#$ -e  /scratch/scratch/ucabtok/130515_nmb_main/000/output/

WDIR=~/Scratch/130515_nmb_main/000/
SCP_OUT=kacprzak@star.ucl.ac.uk:/import/zupcx32/kacprzak/projects/$WDIR/results/
FILENAME_CONFIG=nmb_main.real.test.yaml
DIR_RESULTS=results
DIR_OUTPUT=output
N_OBJ=8
DIR_BIN=~/nmb/nmb_main/
FILENAME_GIT=gitversion.txt

# -------------- do not modify beyond this point ---------------

# calculate object number
TASK_ID=$((SGE_TASK_ID-1))
N_TASKS=$((SGE_TASK_LAST-SGE_TASK_FIRST+1))
OBJ_NUM=$((TASK_ID*N_OBJ))

# make results dicts
mkdir "$WDIR/$DIR_RESULTS"
mkdir "$WDIR/$DIR_OUTPUT"

# get the git version
cd $DIR_BIN
git rev-parse HEAD > "$WDIR/$FILENAME_GIT"

# we always work in tempdir
cd $TMPDIR

# greeting
echo $TASK_ID `date` "hpc is alive"
echo WDIR is $WDIR
echo SGE_TASK_LAST $SGE_TASK_LAST
echo SGE_TASK_FIRST $SGE_TASK_FIRST
echo TASK_ID $TASK_ID
echo JOB_NAME  $JOB_NAME 
echo $OBJ_NUM
echo $N_OBJ

# load modules
echo $TASK_ID `date` "loading modules"
source ~/source_all.sh

# create command
echo $TASK_ID `date` "creating command"
CMD="python $DIR_BIN/nmb_main.py run $WDIR/$FILENAME_CONFIG -v 2 --obj_num $OBJ_NUM --nimages $N_OBJ"
echo $CMD
echo $TASK_ID `date` "running command"
$CMD
echo `date` "copying"
cp $TMPDIR/results.* $WDIR/$DIR_RESULTS/

# Edit and uncomment to copy results to home.
scp $TMPDIR/results.* $SCP_OUT/

echo `date` "exiting submission script"

