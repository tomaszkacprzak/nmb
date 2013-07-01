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

on @legion 
python ~/nmb/nmb_main/nmb_main.py nmb_main.real.optimize.yaml --filepath_truth ~/nmb/nmb_main/truth.optimize.cat -v3

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
m1 = -0.0395     +/-  0.0015
m2 = -0.0396     +/-  0.0014
c1 = -0.0013     +/-  0.0001
c2 = -0.0014     +/-  0.0001

set up 202_bfit_noisy run

from legion, got the commit number for im3shape: 3be9e7bf10aa (levmar_eps5) tip
updating im3shape to this version on precise

test the code with 
precise$ python /home/tomek/Work/code/nmb/nmb_main/nmb_main.py ~/Work/code/nmb/nmb_main/nmb_main.bfit.noisy.yaml -v3 --filepath_truth bfit.nmb_main.real.fits --nimages 64
ran through

lunching on legion
run through

130617 trying to merge results for bfit -- qstat shows only 1/2 done
merge the results in 202_bfit_noisy
zupcx32$ python ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.bfit.noisy.yaml --filepath_truth truth.25880.fits -v3

130618
writing functions to get bias for bins
python ~/code/nmb/nmb_main/nmb_main_plots.py  nmb_main.real.noisy.yaml --filepath_results results.nmb_main.real.noisy.fits --filepath_truth truth.25880.fits -v3 --filepath_acs /home/kacprzak/code/nmb/nmb_main/data/cosmos_acs_shera_may2011.fits.gz

130619
The real noisy results seem to be solid - I checked the errorbars and it looks OK
The results are:

Now run for the noiseless case:

    zupcx32$ python ~/code/nmb/nmb_main/nmb_main_plots.py  nmb_main.real.yaml --filepath_results results.nmb_main.real.pp --filepath_truth truth.26000.pp -v2 --filepath_acs /home/kacprzak/code/nmb/nmb_main/data/cosmos_acs_shera_may2011.fits.gz
    loading truth.26000.pp
    loaded truth.26000.pp correctly, got 1664000 rows
    loading results.nmb_main.real.pp
    loaded results.nmb_main.real.pp correctly, got 1664000 rows
    getting results for all galaxies 1664000
    m1 = -0.0034     +/-  0.0002
    m2 = -0.0038     +/-  0.0002
    c1 =  0.0001     +/-  0.0000
    c2 =  0.0001     +/-  0.0000
    getting results for bins
    redshift bin 0, number of galaxies in sample 435520
    m1 = -0.0049     +/-  0.0002
    m2 = -0.0057     +/-  0.0002
    c1 =  0.0002     +/-  0.0000
    c2 =  0.0002     +/-  0.0000
    redshift bin 1, number of galaxies in sample 506176
    m1 = -0.0039     +/-  0.0002
    m2 = -0.0042     +/-  0.0001
    c1 =  0.0001     +/-  0.0000
    c2 =  0.0001     +/-  0.0000
    redshift bin 2, number of galaxies in sample 376704
    m1 = -0.0027     +/-  0.0002
    m2 = -0.0028     +/-  0.0002
    c1 =  0.0000     +/-  0.0000
    c2 =  0.0001     +/-  0.0000
    redshift bin 3, number of galaxies in sample 303360
    m1 = -0.0013     +/-  0.0002
    m2 = -0.0017     +/-  0.0002
    c1 = -0.0000     +/-  0.0000
    c2 = -0.0000     +/-  0.0000
    redshift bin 4, number of galaxies in sample 40064
    m1 = -0.0014     +/-  0.0009
    m2 = -0.0023     +/-  0.0017
    c1 = -0.0001     +/-  0.0001
    c2 =  0.0002     +/-  0.0001

for the noisy real

    zupcx32$ python ~/code/nmb/nmb_main/nmb_main_plots.py  nmb_main.real.noisy.yaml --filepath_results results.nmb_main.real.noisy.fits --filepath_truth truth.25880.fits -v2 --filepath_acs /home/kacprzak/code/nmb/nmb_main/data/cosmos_acs_shera_may2011.fits.gz
    loading truth.25880.fits
    loaded truth.25880.fits correctly, got 1656320 rows
    loading results.nmb_main.real.noisy.fits
    loaded results.nmb_main.real.noisy.fits correctly, got 1656320 rows
    getting results for all galaxies 1656320
    m1 = -0.0396     +/-  0.0018
    m2 = -0.0397     +/-  0.0016
    c1 = -0.0013     +/-  0.0002
    c2 = -0.0014     +/-  0.0001
    getting results for bins
    redshift bin 0, number of galaxies in sample 431680
    m1 = -0.0209     +/-  0.0027
    m2 = -0.0240     +/-  0.0028
    c1 = -0.0011     +/-  0.0002
    c2 = -0.0012     +/-  0.0002
    redshift bin 1, number of galaxies in sample 501504
    m1 = -0.0424     +/-  0.0027
    m2 = -0.0380     +/-  0.0023
    c1 = -0.0011     +/-  0.0002
    c2 = -0.0014     +/-  0.0002
    redshift bin 2, number of galaxies in sample 373056
    m1 = -0.0548     +/-  0.0035
    m2 = -0.0529     +/-  0.0030
    c1 = -0.0015     +/-  0.0003
    c2 = -0.0015     +/-  0.0003
    redshift bin 3, number of galaxies in sample 300608
    m1 = -0.0392     +/-  0.0036
    m2 = -0.0407     +/-  0.0035
    c1 = -0.0013     +/-  0.0003
    c2 = -0.0017     +/-  0.0003
    redshift bin 4, number of galaxies in sample 39680
    m1 = -0.0577     +/-  0.0098
    m2 = -0.0483     +/-  0.0069
    c1 = -0.0026     +/-  0.0008
    c2 = -0.0010     +/-  0.0006

merged run 202_bfit_noisy

    zupcx32$ python ~/code/nmb/nmb_main/nmb_main_plots.py  nmb_main.bfit.noisy.yaml --filepath_results results.nmb_main.bfit.noisy.fits --filepath_truth truth.25880.fits -v2 --filepath_acs /home/kacprzak/code/nmb/nmb_main/data/cosmos_acs_shera_may2011.fits.gz
    loading truth.25880.fits
    loaded truth.25880.fits correctly, got 1656320 rows
    loading results.nmb_main.bfit.noisy.fits
    loaded results.nmb_main.bfit.noisy.fits correctly, got 1656320 rows
    getting results for all galaxies 1656320
    m1 = -0.0364     +/-  0.0012
    m2 = -0.0367     +/-  0.0009
    c1 = -0.0015     +/-  0.0001
    c2 = -0.0014     +/-  0.0001
    getting results for bins
    redshift bin 0, number of IDs in ACS 8164, number of galaxies in sample 427392
    m1 = -0.0251     +/-  0.0026
    m2 = -0.0168     +/-  0.0027
    c1 = -0.0016     +/-  0.0002
    c2 = -0.0015     +/-  0.0002
    redshift bin 1, number of IDs in ACS 9109, number of galaxies in sample 497792
    m1 = -0.0347     +/-  0.0025
    m2 = -0.0401     +/-  0.0028
    c1 = -0.0011     +/-  0.0002
    c2 = -0.0016     +/-  0.0002
    redshift bin 2, number of IDs in ACS 6725, number of galaxies in sample 370624
    m1 = -0.0482     +/-  0.0036
    m2 = -0.0451     +/-  0.0037
    c1 = -0.0015     +/-  0.0003
    c2 = -0.0011     +/-  0.0003
    redshift bin 3, number of IDs in ACS 5435, number of galaxies in sample 298496
    m1 = -0.0362     +/-  0.0036
    m2 = -0.0412     +/-  0.0043
    c1 = -0.0020     +/-  0.0003
    c2 = -0.0015     +/-  0.0004
    redshift bin 4, number of IDs in ACS 739, number of galaxies in sample 39360
    m1 = -0.0407     +/-  0.0163
    m2 = -0.0604     +/-  0.0116
    c1 = -0.0031     +/-  0.0014
    c2 = -0.0010     +/-  0.0010

these were done using commit 0eec29501273fdedd6691a1d66b8882e1d3437bc

now rewriting binning functions to be more general

    zupcx32$ python ~/code/nmb/nmb_main/nmb_main_plots.py  nmb_main.real.noisy.yaml --filepath_results results.nmb_main.real.noisy.fits --filepath_truth truth.25880.fits -v2 --filepath_acs /home/kacprzak/code/nmb/nmb_main/data/cosmos_acs_shera_may2011.fits.gz --filepath_stats /home/kacprzak/code/nmb/nmb_main/data/stats.nmb_main.real.pp

the functions now use fixed results filenames, call
    
    precise$ python ~/Work/code/nmb/nmb_main/nmb_main_plots.py -v 2

commit 9c0c71161ba00d982b15e6a593ae9cbe78745ce8
on zupcx32 using python 2.7    

    zupcx32$ python $CODE/nmb/nmb_main/nmb_main_plots.py -v 2

commit 956222c 

130621 found the indexing bug and made a lot of new plots    

    precise$ python ~/Work/code/nmb/nmb_main/nmb_main_plots.py -v 2

    submitting qsub nmb_main.legion.bfit.noisy.run2.sh

130622 using commit e8509a03ad6bab68f5813967fc313068c604f241

    now the legion run is finished --> incorporating the new results


    Turns out that we are missing around 800 files, which makes 8000 galaxies -- really bad. Remember to increase the clock time to 6hrs!

    python ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.bfit.noisy.yaml --filepath_truth truth.25880.fits -v3 --n_reps 2

    results  n    3312640 first 1000270000 last 0
    truth    n    1656320 first 1000270000 last 742430707
    getting number of matches ...
    number of matches 2800640
    number of matches 3312640
    saving tables ...
    Overwriting existing file 'results.nmb_main.bfit.noisy.fits'.
    table saved results.nmb_main.bfit.noisy.fits correctly, got 1 rows
    results saved results.nmb_main.bfit.noisy.fits correctly, got 3312640 rows
    truth   n    3312640 first 1000270000 last 742430707
    results n    1656320 first 1000270000 last 742430707

    looks like everything is matched correctly!

    using 
    testSaveBiasForBins()

using the old results file filepath_results_bfit_noisy = 'results.nmb_main.bfit.noisy.fits', bin 2:

m1 =  0.0357     +/-  0.0033
m2 =  0.0348     +/-  0.0033
c1 = -0.0018     +/-  0.0003
c2 = -0.0017     +/-  0.0003

using the new results file filepath_results_bfit_noisy = 'results.nmb_main.bfit.noisy.rep2.fits', bin 2:

almost the same with +/-  0.0025 errors


made new plots for m1 m2 vs model bias -- need more statistics

second real_noisy is finished

    python ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.bfit.noisy.yaml --filepath_truth truth.25880.fits -v3 --n_reps 2

130107 got so far so many results:
bfit noisy 4416 - 2 reps
real noisy 5839 - 3 reps

merging 302 results

        python ~/code/nmb/nmb_main/nmb_main_analyse.py mergeResults nmb_main.real.noisy.yaml --filepath_truth truth.25880.fits -v3 --n_reps 3

merged fine:

    7763 file results/results.nmb_main.real.noisy.000004968320.cat n_gals 640
    results  n    4968960 first 1000270000 last 742430707
    truth    n    1656320 first 1000270000 last 742430707
    getting number of matches ...
    number of matches 3722240
    number of matches 4968960
    saving tables ...
    Overwriting existing file 'results.nmb_main.real.noisy.fits'.
    table saved results.nmb_main.real.noisy.fits correctly, got 1 rows
    results saved results.nmb_main.real.noisy.fits correctly, got 4968960 rows
    results   n    4968960 first 1000270000 last 742430707
    truth     n    1656320 first 1000270000 last 742430707

plots:

    zupcx32$ python $CODE/nmb/nmb_main/nmb_main_plots.py plotBiasForBins -v 2

    looks OK, but some backend problems, I will try on laptop

in the meanwhile submitted the third bfit run:
[ucabtok@login01 202_bfit_noisy]$ qsub nmb_main.legion.bfit.noisy.run3.sh



