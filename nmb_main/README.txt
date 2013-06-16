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
zupcx32$ python ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.real.noisy.yaml --filepath_truth truth.25880.fits -v3
table saved results.nmb_main.real.noisy.fits correctly, got 1 rows
results saved results.nmb_main.real.noisy.fits correctly, got 1656320 rows
truth   n    1656320 first 1000270000 last 742430707
results n    1656320 first 1000270000 last 742430707

zupcx32$ python ~/code/nmb/nmb_main/nmb_main_analyse.py getTotalBias nmb_main.real.noisy.yaml --filepath_truth truth.25880.fits -v3
got very bad m results:
m1 = -0.2797     +/-  0.0012
m2 = -0.2814     +/-  0.0012
c1 =  0.0006     +/-  0.0001
c2 =  0.0003     +/-  0.0001
will create a test scrit to tweak levmar settings

cd @zupcx32 130515_nmb_main/002/
python ~/code/nmb/nmb_main/nmb_main.py ~/code/nmb/nmb_main/nmb_main.real.optimize.yaml --filepath_truth ~/code/nmb/nmb_main/truth.optimize.cat

changed the system so that the images are normalised to unit flux (and noise sigma = 1), so that the optimization is more stable
launched a run 302_real_noisy, with results from the optimization

130616 can't believe that it has been a almost 3 weeks since I have been working on this!

merge the results in 302_real_noisy
zupcx32$ python ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.real.noisy.yaml --filepath_truth truth.25880.fits -v3
...
saving tables ...
table saved results.nmb_main.real.noisy.fits correctly, got 1 rows
results saved results.nmb_main.real.noisy.fits correctly, got 1656320 rows
truth   n    1656320 first 1000270000 last 742430707
results n    1656320 first 1000270000 last 742430707

get total bias
zupcx32$ python ~/code/nmb/nmb_main/nmb_main_analyse.py getTotalBias nmb_main.real.noisy.yaml --filepath_truth truth.25880.fits -v3
woah it actually worked and gives 'sensible' results!
m1 = -0.0395 	 +/-  0.0015
m2 = -0.0396 	 +/-  0.0014
c1 = -0.0013 	 +/-  0.0001
c2 = -0.0014 	 +/-  0.0001

set up 202_bfit_noisy run

test the code with 
zupcx32$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.bfit.noisy.yaml -v3 --filepath_truth bfit.nmb_main.real.fits --nimages 64



