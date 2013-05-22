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
import scipy.io
import cPickle as pickle

# global results_array

dtype_table_truth   = { 'names'  : ['id_unique','id_cosmos','g1','g2','angle','id_angle','id_shear' , 'zphot'],
                        'formats': ['i8']*2 + ['f4']*3 + ['i4']*2 + ['f4']*1 }

dtype_table_results = { 'names'   : ['identifier','likelihood','time_taken','x0','y0','e1','e2','radius','fwhm','bulge_flux','disc_flux','flux_ratio','signal_to_noise','min_residuals','max_residuals','model_min','model_max','number_of_likelihood_evals','number_of_iterations','reason_of_termination'],
                        'formats' : ['i4'] + ['f4']*16 + ['i4']*3 }            

dtype_table_stats =  { 'names'   : ['index', 'cosmos_id','m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std' , 'hlr' , 'rgp' , 'snr'],
                        'formats' : ['i4']*2 + ['f4']*11 }            

filename_results_pickle = 'results.nmb_main.real.pp'
filename_stats_pickle = 'stats.nmb_main.real.pp'


NO_RESULT_FLAG = 666

def _getLineFit(x,y,sig):
        """
        @brief get linear least squares fit with uncertainity estimates
        y(X) = b*X + a
        see numerical recipies 15.2.9
        @param X    function arguments 
        @param y    function values
        @param sig  function values standard deviations   
        @return a - additive 
        @return b - multiplicative
        @return C - covariance matrix
        """
        
        import numpy

        invsig2 = sig**-2;
        
        S  = numpy.sum(invsig2)
        Sx = numpy.inner(x,invsig2)
        Sy = numpy.inner(y,invsig2)
        Sxx = numpy.inner(invsig2*x,x)
        Sxy = numpy.inner(invsig2*x,y)
        
        D = S*Sxx - Sx**2
        a = (Sxx*Sy - Sx*Sxy)/D
        b = (S*Sxy  - Sx*Sy)/D
        
        Cab = numpy.zeros((2,2))
        Cab[0,0] = Sxx/D
        Cab[1,1] = S/D
        Cab[1,0] = -Sx/D
        Cab[0,1] = Cab[1,0]
        
        return a,b,Cab

def loadResultsArray():
    """
    Loads the resutls array from file global filename_results_pickle, unless there is an ipython global with results_array.
    @return results_array
    """

    # this loads the 
    global results_array
    
    file_pickle = open(filename_results_pickle)
    try:
        results_array
    except:
        logger.info('loading %s' % filename_results_pickle)
        results_array = pickle.load(file_pickle)
        file_pickle.close()
    else:
        logger.info('using preloaded results_array')
    
    logger.info('loaded %s correctly, got %d rows' % (filename_results_pickle,results_array.shape[0]))

    return results_array

def loadTruthArray():


    # this loads the 
    global truth_array

    filename_pickle = config['args'].filepath_truth
    file_pickle = open(filename_pickle)
 
    try:
        truth_array
    except:
        logger.info('loading %s' % config['args'].filepath_truth)
        truth_array = pickle.load(file_pickle)
    else:
        logger.info('using preloaded truth_array')
    
    logger.info('loaded %s correctly, got %d rows' % (config['args'].filepath_truth,truth_array.shape[0]))

    return truth_array

def getBiasForEachGal():

    # get the results and truth
    results_array = loadResultsArray()
    truth_array   = loadTruthArray()

    cosmos_ids = set(truth_array['id_cosmos'])
    n_per_cosmos = sum(truth_array['id_cosmos'] == truth_array[0]['id_cosmos'])
    print "n_per_cosmos" , n_per_cosmos
    n_shears = len(set(truth_array['id_shear']))
    n_angles = len(set(truth_array['id_angle']))

    logger.info('using set of %d cosmos galaxies ' % len(cosmos_ids) )

    results_stats = numpy.zeros(1,dtype=dtype_table_stats)

    n_gals_lost = 0

    for i,cid in enumerate(cosmos_ids):

        shears_true_g1 = []
        shears_true_g2 = []
        shears_mean_g1 = []
        shears_mean_g2 = []
        shears_stdm_g1 = []
        shears_stdm_g2 = []
        rgp_list       = []
        snr_list       = []
        hlr_list       = []

        for gid in range(n_shears):

            n_start = i*n_per_cosmos + gid     * n_angles
            n_end   = i*n_per_cosmos + (gid+1) * n_angles

            truth_current   = truth_array[n_start:n_end]
            results_current = results_array[n_start:n_end]

            if any(results_current['e1'] == NO_RESULT_FLAG):
                logger.error('% 8d %10d -- missing angles for shear %d' % (i,cid,gid))
                continue

            mean_g1 = numpy.mean(results_current['e1'])      
            mean_g2 = numpy.mean(results_current['e2'])      
            
            if numpy.isnan(mean_g1) or numpy.isnan(mean_g2):
                logger.error('% 8d %10d -- shear is nan for %d' % (i,cid,gid))
                continue

            shears_true_g1.append(      truth_current['g1'][0]                              )
            shears_true_g2.append(      truth_current['g2'][0]                              )
            shears_mean_g1.append(      mean_g1                                             )
            shears_mean_g2.append(      mean_g2                                             )
            shears_stdm_g1.append(      numpy.std(results_current['e1'],ddof=1)             )
            shears_stdm_g2.append(      numpy.std(results_current['e2'],ddof=1)             )
            rgp_list.append(            numpy.mean(results_current['fwhm'])                 )  
            snr_list.append(            numpy.mean(results_current['signal_to_noise'])      )  
            hlr_list.append(            numpy.mean(results_current['radius'])               )  


        n_valid_shears = len(shears_mean_g1)
        if n_valid_shears < 5:
            n_gals_lost+=1
            logger.error('% 8d %10d -- not enough shears for galaxy : %d , so far lost %d' % (i,cid,n_valid_shears,n_gals_lost))
            continue
        else:
            bias_g1 = numpy.array(shears_mean_g1,dtype=numpy.float64) - numpy.array(shears_true_g1,dtype=numpy.float64)
            bias_g2 = numpy.array(shears_mean_g2,dtype=numpy.float64) - numpy.array(shears_true_g2,dtype=numpy.float64)
            g1_true = numpy.array(shears_true_g1,dtype=numpy.float64) 
            g2_true = numpy.array(shears_true_g2,dtype=numpy.float64)
            bias_g1_stdm = numpy.ones(bias_g1.shape) 
            bias_g2_stdm = numpy.ones(bias_g2.shape)

            c1,m1,cov1 = _getLineFit(g1_true,bias_g1,bias_g1_stdm)
            c2,m2,cov2 = _getLineFit(g2_true,bias_g2,bias_g2_stdm)
            m1_std = numpy.sqrt(cov1[1,1])
            m2_std = numpy.sqrt(cov2[1,1])
            c1_std = numpy.sqrt(cov1[0,0])
            c2_std = numpy.sqrt(cov2[0,0])
    
            rgp = numpy.mean(rgp_list)
            snr = numpy.mean(snr_list)
            hlr = numpy.mean(hlr_list)

            if any(numpy.isnan([m1, m2, c1, c2])):
                n_gals_lost+=1
                logger.error('% 8d %10d -- nans in bias line fit params : %d , so far lost %d' % (i,cid,n_valid_shears,n_gals_lost))
                continue

            # dtype_table_stats =  { 'names'   : ['index', 'cosmos_id','m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std' , 'hlr' , 'rgp' , 'snr'],
            results_stats_row = numpy.array([(i,cid,m1,m2,m1_std,m2_std,c1,c2,c1_std,c2_std,hlr,rgp,snr)],dtype=dtype_table_stats)
            results_stats = numpy.concatenate((results_stats,results_stats_row))


        logger.debug('%8d galaxy %8d valid shears %d, so far got %8d valid gals' % (i,cid,n_valid_shears,results_stats.shape[0]))

    # remove the first zero row
    results_stats=results_stats[1:]
    n_stats = len(results_stats)
    logger.info('results_stats has %d rows' % n_stats)

    file_pickle = open(filename_stats_pickle,'w')
    pickle.dump(results_stats,file_pickle,protocol=2)
    file_pickle.close()
    logger.info('saved %s' % filename_stats_pickle) 
    
def mergeResutls():

    truth_array   = loadTruthArray()

    # get the wildcard for the files - join all files in a big array
    filecard_resutls = os.path.join(config['args'].dirpath_results,'results.nmb_main.real.*.cat')
    logger.info(filecard_resutls)
    # count files and columns
    
    n_gals_total = config['settings']['n_images']
    n_gals_per_file = 640
    files = [os.path.join(config['args'].dirpath_results,'results.nmb_main.real.%05d.cat' % n) for n in range(0,n_gals_total,n_gals_per_file)]
    n_files = len(files)
    n_meas = len(dtype_table_results['names'])
    logger.info('got %d file names' % n_files)

    # initialise the empty array
    results_array = numpy.zeros(n_gals_total,dtype=dtype_table_results)
    results_error = numpy.zeros(1,dtype=dtype_table_results)
    results_error['e1'] = NO_RESULT_FLAG 
    results_error['e2'] = NO_RESULT_FLAG
    
    # truth index
    ti = 0

    # loop over files
    for fi,file_results in enumerate(files):

        # load file
        if os.path.isfile(file_results):

            results = numpy.loadtxt(file_results,dtype=dtype_table_results)
            n_gals_in_file = len(results)

            for gi in range(n_gals_per_file):

                if gi < n_gals_in_file:

                    if results[gi]['identifier'] == truth_array[ti]['id_unique']:
                        results_array[ti] = results[gi]
                    else:
                        logger.error('gi in bounds, but no entry for file %s galaxy %d id_result = %d id_unique = %d' % (file_results,gi,results['identifier'][gi],truth_array[ti]['id_unique']))
                        results_array[ti] = results_error
                        retults_array[ti]['identifier'] = truth_array[ti]['id_unique']
                else:
                    logger.error('no entry for file %s galaxy %d ti %d id_unique = %d' % (file_results,gi,ti,truth_array[ti]['id_unique']))
                    results_array[ti] = results_error
                    results_array[ti]['identifier'] = truth_array[ti]['id_unique']

                ti+=1

            logger.info('% 4d\t% 4d\t%50s\t%d galaxies' % (fi,ti,file_results,n_gals_in_file)) 

        else:
            logger.error('file %s not found' % fi)
            for gi in range(n_gals_per_file):
                results_array[ti] = results_error
                ti+=1
     

    # save file pickle - use global filename_pickle
    file_pickle = open(filename_results_pickle,'w')
    logger.info('saving')
    pickle.dump(results_array,file_pickle,protocol=2)
    file_pickle.close()
    # check if OK
    results_array = pickle.load(open(filename_results_pickle))
    logger.info('saved %s correctly, got %d rows' % (filename_results_pickle,results_array.shape[0]))
    file_pickle.close()

def plotBiasHistogram():

    def _binCenters(b):
          c1 = [ b[i-1] + (b[i] - b[i-1])/2. for i in range(1,len(b)) ] 
          return c1


    file_pickle = open(filename_stats_pickle)
    results_stats = pickle.load(file_pickle)
    logger.info('opened stats file %s with %d rows' % (filename_stats_pickle,results_stats.shape[0]))

    n_bins = 32
    m_cut = 0.1
    select1 = numpy.logical_and(results_stats['m1'] < m_cut , results_stats['m1'] > -m_cut)
    select2 = numpy.logical_and(results_stats['m2'] < m_cut , results_stats['m2'] > -m_cut) 
    select = numpy.logical_and(select1,select2)
    results_m1 = results_stats[select]['m1']
    results_m2 = results_stats[select]['m2']
    results_rgp = results_stats[select]['rgp']
    results_snr = results_stats[select]['snr']

    pylab.ion()
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
    pylab.show()
    filename_fig = 'figures/figure.hist.m1m2.real.png'
    pylab.savefig(filename_fig)

    pylab.figure()
    bins_snr = numpy.logspace(2,4,n_bins)
    h,b,p = pylab.hist(results_snr,bins=bins_snr)
    pylab.xlabel('snr')
    pylab.xscale('log')
    pylab.show()    
    filename_fig = 'figures/figure.hist.snr.real.png'
    pylab.savefig(filename_fig)
    
    pylab.figure()

    psf_size =  0.7695
    results_size = results_rgp/psf_size
    bins_size = numpy.linspace(1,4,n_bins)
    h,b,p = pylab.hist(results_size,bins=bins_size)
    pylab.xlabel('rgp/rp')
    pylab.show()    
    filename_fig = 'figures/figure.hist.size.real.png'
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
    filename_fig = 'figures/figure.hist.m1m2vsSNR.real.png'
    pylab.savefig(filename_fig)

    # m1 vs size


    digitized = numpy.digitize(results_size, bins_size)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(0, len(bins_size))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(0, len(bins_size))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(0, len(bins_size))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(0, len(bins_size))]
    
    pylab.figure()
    pylab.errorbar(bins_size,bin_mean1,yerr=bin_stdv2)
    pylab.errorbar(bins_size,bin_mean2,yerr=bin_stdv2)
    pylab.xlabel('size')
    pylab.ylabel('m1')

    filename_fig = 'figures/figure.hist.m1m2vsSIZE.real.png'
    pylab.savefig(filename_fig)

    print len(results_size[results_size > 1.2])


def main():

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    parser.add_argument('filepath_truth', type=str, default=None, help='truth file for the run, overrides the config file (by default is taken from yaml file)')
    parser.add_argument('--dirpath_results', type=str, default='results', help='where the results files are')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')

    args = parser.parse_args()
    args.name_config = os.path.basename(args.filepath_config).replace('.yaml','')

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    global logger , config 
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("nmb_main") 
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