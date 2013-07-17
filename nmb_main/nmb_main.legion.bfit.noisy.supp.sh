#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=6:0:0
#$ -l mem=2G
#$ -t 12944-15532
#$ -N nmb_bfit_noisy
#$ -P CosmicShear
#$ -wd /scratch/scratch/ucabtok/130515_nmb_main/202_bfit_noisy/
#$ -o  /scratch/scratch/ucabtok/130515_nmb_main/202_bfit_noisy/output/
#$ -e  /scratch/scratch/ucabtok/130515_nmb_main/202_bfit_noisy/output/
#Local2Scratch

WDIR=~/Scratch/130515_nmb_main/202_bfit_noisy/
SCP_OUT=kacprzak@star.ucl.ac.uk:/import/zupcx32/kacprzak/projects/130515_nmb_main/202_bfit_noisy/results/
FILENAME_CONFIG=nmb_main.bfit.noisy.yaml
FILENAME_TRUTH=bfit.nmb_main.real.fits
FILENAME_INI=nmb.ini
FILENAME_GIT=gitversion.txt
FILENAME_SUPP=suppl.cat
DIR_RESULTS=results
DIR_OUTPUT=output
N_OBJ=640
DIR_BIN=~/nmb/nmb_main/

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
echo OBJ_NUM $OBJ_NUM
echo N_OBJ $N_OBJ
echo TMPDIR $TMPDIR

# check if files are available
ls $DIR_BIN/nmb_main.py 
ls $WDIR/$FILENAME_CONFIG
ls $WDIR/$FILENAME_INI
ls $WDIR/$FILENAME_TRUTH
ls $WDIR/$FILENAME_SUPP

# load modules
echo $TASK_ID `date` "loading modules"
source ~/source_all.sh
source ~/source_paths.sh

# create command
echo $TASK_ID `date` "creating command"
CMD="python $DIR_BIN/nmb_main.py $WDIR/$FILENAME_CONFIG --filepath_ini $WDIR/$FILENAME_INI -v 1 --obj_num $OBJ_NUM --nimages $N_OBJ --filepath_truth $WDIR/$FILENAME_TRUTH --filepath_idslist $WDIR/$FILENAME_SUPP"
echo $CMD
echo $TASK_ID `date` "running command"
$CMD
echo `date` "copying"
cp $TMPDIR/results.* $WDIR/$DIR_RESULTS/

# Edit and uncomment to copy results to home.
scp $TMPDIR/results.* $SCP_OUT/

echo `date` "exiting submission script"
