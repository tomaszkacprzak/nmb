# This file is a log of the nmb experiment. Hopefully it will allow to reproduce the procedure of getting the results for the paper.

# 130516 create the truth tables
# precise$ python nmb_truth.py --filepath_ids galids.2.cat     --filepath_out truth.test.cat  --n_angles 4 --shears 0.05 
# precise$ python nmb_truth.py --filepath_ids galids.26000.cat --filepath_out truth.26000.cat --n_angles 8 --shears 0.1 

# 130517 run the test scripts
# python nmb_main.py run nmb_main.real.test.yaml -v 2
# also runs on legion
# nmb_main.legion.test.sh

# 130515 check if the ini settings are reproducing galsim images well  
# run compare_im3shape_galsim.py

# 130515 run on legion

# 130517 got results for 101

# 130518 merge resutls
# mergeResutls

# 130521 produce first plots
# ipython run -i /home/tomek/Work/code/nmb/nmb_main/nmb_main_plots.py nmb_main.real.yaml truth.26000.pp -v 2

# 130522 update the table with redshifts:
# precise$ python nmb_truth.py --filepath_ids galids.26000.cat  --filepath_z cosmos_acs_shera_may2011.fits.gz --filepath_out truth.26000.cat  --n_angles 8 --shears 0.1 

# 130523 
# merge
# run ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.real.yaml --filepath_truth truth.26000.pp -v3
# got stats
# run ~/code/nmb/nmb_main/nmb_main_analyse.py getBiasForEachGal nmb_main.real.yaml -v2

# 130524
# create truth table for the bfit table results.bfit.nmb_main.real.fits  truth.bfit.26000.fits in 150515_nmb_main/101
# copy the above to 150515_nmb_main/003

# produce plots with hlr
# In [23]: run -i /home/tomek/Work/code/nmb/nmb_main/nmb_main_plots.py nmb_main.real.yaml  -v 2

# 130525
# updated im3shape to revision changeset: 177:608e4414bffe
# run two test scripts:
# precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.real.test.yaml -v3 --filepath_truth ~/Work/code/nmb/nmb_main/truth.26000.cat 
# precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.bfit.test.yaml -v3 
# for the bfit, I obtained model bias of order 0.0001

# created data directory and moved there all the big files

# 130526 
# @precise 150515_nmb_main/003
# precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.bfit.test.yaml -v3 --filepath_truth ~/Work/code/nmb/nmb_main/data/bfit.nmb_main.real.fits -snr 20.0
# precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.real.test.yaml -v3 --filepath_truth ~/Work/code/nmb/nmb_main/data/truth. -snr 20.0

# get the bfit catalogs
# precise$ python ~/Work/code/nmb/nmb_main/nmb_main_analyse.py createBFITsample nmb_main.real.yaml 
# the file is truth.25880.fits

# test the noisy scripts
# precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.bfit.noisy.yaml -v3 --filepath_truth ~/Work/code/nmb/nmb_main/data/bfit.nmb_main.real.fits --nimages 64
# precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.real.noisy.yaml -v3 --filepath_truth ~/Work/code/nmb/nmb_main/data/truth.25880.fits --nimages 64

# test @legion
# time python ~/nmb/nmb_main/nmb_main.py ~/nmb/nmb_main/nmb_main.bfit.noisy.yaml --filepath_truth ~/nmb/nmb_main/data/bfit.nmb_main.real.fits --nimages 64 --filepath_ini ~/nmb/nmb_main/nmb.ini
# time python ~/nmb/nmb_main/nmb_main.py ~/nmb/nmb_main/nmb_main.real.noisy.yaml --filepath_truth ~/nmb/nmb_main/data/truth.25880.fits --nimages 64 --filepath_ini ~/nmb/nmb_main/nmb.ini

# qrsh -l mem=2G,h_rt=1:0:0 -P CosmicShear
# time @legion 120s / 64 gals = 1200s / 640 gals = 20 min per job -- set the wallclock to 2h

# create 130515_nmb_main/201 -- real
# create 130515_nmb_main/301 -- bfit

# submit legion 301 job -- got it back

# 130530
# run ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.real.noisy.yaml --filepath_truth truth.25880.fits -v3


