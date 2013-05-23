import os
import pyfits
import logging
import sys
import numpy
import sys
import math
import argparse
import yaml
import galsim
import copy
import datetime  
import pylab
import glob
import cPickle as pickle

# global results_array

dtype_table_truth   = { 'names'  : ['id_unique','id_cosmos','g1','g2','angle','id_angle','id_shear' , 'zphot'],
                        'formats': ['i8']*2 + ['f4']*3 + ['i4']*2 + ['f4']*1 }

dtype_table_results = { 'names'   : ['identifier','likelihood','time_taken','x0','y0','e1','e2','radius','fwhm','bulge_flux','disc_flux','flux_ratio','signal_to_noise','min_residuals','max_residuals','model_min','model_max','number_of_likelihood_evals','number_of_iterations','reason_of_termination'],
                        'formats' : ['i8'] + ['f4']*16 + ['i4']*3 }           

dtype_table_stats =  { 'names'   : ['index', 'cosmos_id' , 'zphot' ,'m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std' , 'hlr' , 'rgp' , 'snr'],
                        'formats' : ['i8']*2 + ['f4']*12 }             

NO_RESULT_FLAG = 666

def plotBiasHistogram():

    def _binCenters(b):
          c1 = [ b[i-1] + (b[i] - b[i-1])/2. for i in range(1,len(b)) ] 
          return c1


    file_pickle = open(args.filepath_stats)
    results_stats = pickle.load(file_pickle)
    logger.info('opened stats file %s with %d rows' % (args.filepath_stats,results_stats.shape[0]))

    n_bins = 32
    m_cut = 0.1
    select1 = numpy.logical_and(results_stats['m1'] < m_cut , results_stats['m1'] > -m_cut)
    select2 = numpy.logical_and(results_stats['m2'] < m_cut , results_stats['m2'] > -m_cut) 
    select = numpy.logical_and(select1,select2)
    select = numpy.logical_and(select , results_stats['zphot'] > 0) 
    results_m1 = results_stats[select]['m1']
    results_m2 = results_stats[select]['m2']
    results_rgp = results_stats[select]['rgp']
    results_snr = results_stats[select]['snr']
    results_zphot = results_stats[select]['zphot']

    # pylab.ion()
    pylab.close('all')

    # m1 m2 hist    

    pylab.figure()
    pylab.clf()
    h1,b1 = pylab.histogram(results_m1,n_bins)
    h2,b2 = pylab.histogram(results_m2,n_bins)
    pylab.plot(_binCenters(b1),h1,'r+-')
    pylab.plot(_binCenters(b2),h2,'bx-')
    pylab.xlabel('m1 , m2')
    print "mean" , numpy.mean(results_m1)
    print "mean" , numpy.mean(results_m2)
    print "std" , numpy.std(results_m1,ddof=1)
    print "std" , numpy.std(results_m2,ddof=1)
    print 'nans?' , any(numpy.isnan(results_m1))
    # pylab.show()
    filename_fig = 'figures/figure.hist.m1m2.real.png'
    pylab.savefig(filename_fig)

    pylab.figure()
    bins_snr = numpy.logspace(2,4,n_bins/2)
    h,b,p = pylab.hist(results_snr,bins=bins_snr)
    pylab.xlabel('snr')
    pylab.xscale('log')
    # pylab.show()    
    filename_fig = 'figures/figure.hist.snr.real.png'
    pylab.savefig(filename_fig)
    
    pylab.figure()
    psf_size =  0.7695
    results_size = results_rgp/psf_size
    bins_size = numpy.linspace(1,4,n_bins)
    h,b,p = pylab.hist(results_size,bins=bins_size)
    pylab.xlabel('rgp/rp')
    # pylab.show()    
    filename_fig = 'figures/figure.hist.size.real.png'
    pylab.savefig(filename_fig)

    # redshift distribution

    pylab.figure()
    bins_redshift = [0 , 0.35, 0.6, 0.8, 1.1 , 1.5]
    h,b,p = pylab.hist(results_zphot,bins=n_bins)
    # h,b,p = pylab.hist(results_zphot,bins=bins_redshift)
    pylab.xlabel('photo-z')
    # pylab.show()    
    filename_fig = 'figures/figure.hist.redshift.real.png'
    pylab.savefig(filename_fig)





    # m1 vs snr

    digitized = numpy.digitize(results_snr, bins_snr)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(0, len(bins_snr))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(0, len(bins_snr))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(0, len(bins_snr))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(0, len(bins_snr))]

    
    pylab.figure()
    pylab.errorbar(bins_snr,bin_mean1,yerr=bin_stdv1)
    pylab.errorbar(bins_snr,bin_mean2,yerr=bin_stdv2)
    pylab.xscale('log')
    pylab.xlabel('snr')
    pylab.ylabel('m1')
    filename_fig = 'figures/figure.hist.M1M2vsSNR.real.png'
    pylab.savefig(filename_fig)

    # m1 vs size


    digitized = numpy.digitize(results_size, bins_size)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(0, len(bins_size))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(0, len(bins_size))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(0, len(bins_size))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(0, len(bins_size))]
    
    pylab.figure()
    pylab.errorbar(bins_size,bin_mean1,yerr=bin_stdv1)
    pylab.errorbar(bins_size,bin_mean2,yerr=bin_stdv2)
    pylab.xlabel('size')
    pylab.ylabel('m1')

    filename_fig = 'figures/figure.hist.M1M2vsSIZE.real.png'
    pylab.savefig(filename_fig)

    # snr vs redshifts

    pylab.figure() 
    digitized = numpy.digitize(results_zphot, bins_redshift)
    for i in range(1, len(bins_redshift)):
        h,b = pylab.histogram(results_snr[digitized == i],bins=bins_snr,normed=True)
        pylab.plot(h,'x-',label='redshift bin %2.2f' % bins_redshift[i])

    pylab.legend()
    pylab.xlabel('snr')
    pylab.ylabel('normalised hist')




    filename_fig = 'figures/figure.hist.SNRvsREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # size vs redshifts

    pylab.figure() 
    digitized = numpy.digitize(results_zphot, bins_redshift)
    for i in range(1, len(bins_redshift)):
        h,b = pylab.histogram(results_size[digitized == i],bins=bins_size,normed=True)
        pylab.plot(h,'x-',label='redshift bin %2.2f' % bins_redshift[i])

    pylab.legend()
    pylab.xlabel('Rgp/Rp')
    pylab.ylabel('normalised hist')

    filename_fig = 'figures/figure.hist.SIZEvsREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # redshifts vs redshifts

    # pylab.figure() 
    # digitized = numpy.digitize(results_zphot, bins_redshift)
    # for i in range(1, len(bins_redshift)):
    #     h,b = pylab.histogram(results_zphot[digitized == i],bins=bins_redshift,normed=True)
    #     pylab.plot(h,'x-',label='redshift bin %2.2f' % bins_redshift[i])

    # pylab.legend()
    # pylab.xlabel('zphot - control')
    # pylab.ylabel('normalised hist')

    # m1 vs redshifts

    digitized = numpy.digitize(results_zphot, bins_redshift)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(1, len(bins_redshift))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(1, len(bins_redshift))]
    bin_stdm1 = [results_m1[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m1[digitized == i])) for i in range(1, len(bins_redshift))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(1, len(bins_redshift))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(1, len(bins_redshift))]
    bin_stdm2 = [results_m2[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m2[digitized == i])) for i in range(1, len(bins_redshift))]
    
    pylab.figure()
    pylab.errorbar(bins_redshift[1:],bin_mean1,yerr=bin_stdv1)
    pylab.errorbar(bins_redshift[1:],bin_mean2,yerr=bin_stdv2)
    pylab.errorbar(bins_redshift[1:],bin_mean1,yerr=bin_stdm1)
    pylab.errorbar(bins_redshift[1:],bin_mean2,yerr=bin_stdm2)
    pylab.xlabel('redshift')
    pylab.ylabel('m1')
    
    filename_fig = 'figures/figure.hist.M1M2vsREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    print len(results_stats)
    print len(results_size[results_size > 1.2])


def main():

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    parser.add_argument('--filepath_truth', type=str, default='truth.26000.sorted.pp', help='truth file for the run, overrides the config file (by default is taken from yaml file)')
    parser.add_argument('--filepath_stats', type=str, default='stats.nmb_main.real.pp', help='stats file')
    parser.add_argument('--filepath_results', type=str, default='results.nmb_main.real.pp', help='results file')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')

    args = parser.parse_args()
    args.name_config = os.path.basename(args.filepath_config).replace('.yaml','')

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    global logger , config , args
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("nmb_main_plots.py") 
    logger.setLevel(logging_level)
    
    # load the configuration file
    config = yaml.load(open(args.filepath_config,'r')) 
    # store the args in config so it's easier to use them
    config['args'] = args    

    # mergeResutls()
    # getBiasForEachGal()
    plotBiasHistogram()

if __name__ == "__main__":

    main()