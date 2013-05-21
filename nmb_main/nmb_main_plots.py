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

dtype_table_truth   = { 'names'  : ['id_unique','id_cosmos','g1','g2','angle','id_angle','id_shear'],
                        'formats': ['i4']*2 + ['f4']*3 + ['i4']*3 }

dtype_table_results = { 'names'   : ['identifier','likelihood','time_taken','x0','y0','e1','e2','radius','fwhm','bulge_flux','disc_flux','flux_ratio','signal_to_noise','min_residuals','max_residuals','model_min','model_max','number_of_likelihood_evals','number_of_iterations','reason_of_termination'],
                        'formats' : ['i4'] + ['f4']*16 + ['i4']*3 }            

dtype_table_stats =  { 'names'   : ['unique_id','gal_num','m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std'],
                        'formats' : ['i4']*2 + ['f4']*8 }            

filename_pickle = 'results.nmb_main.real.pp'

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
    Loads the resutls array from file global filename_pickle, unless there is an ipython global with results_array.
    @return results_array
    """

    # this loads the 
    global results_array
    
    file_pickle = open(filename_pickle)
    try:
        results_array
    except:
        logger.info('loading %s' % filename_pickle)
        results_array = pickle.load(file_pickle)
        file_pickle.close()
    else:
        logger.info('using preloaded results_array')
    
    logger.info('loaded %s correctly, got %d rows' % (filename_pickle,results_array.shape[0]))

    return results_array

def loadTruthArray():


    # this loads the 
    global truth_array
    
    file_pickle = open(config['args'].filepath_truth)
    try:
        truth_array
    except:
        logger.info('loading %s' % config['args'].filepath_truth)
        truth_array = numpy.loadtxt(config['args'].filepath_truth,dtype=dtype_table_truth)
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

    for i,cid in enumerate(cosmos_ids):

        shears_true_g1 = []
        shears_true_g2 = []
        shears_mean_g1 = []
        shears_mean_g2 = []
        shears_stdm_g1 = []
        shears_stdm_g2 = []

        for gid in range(n_shears):

            n_start = i*n_per_cosmos + gid     * n_angles
            n_end   = i*n_per_cosmos + (gid+1) * n_angles

            truth_current   = truth_array[n_start:n_end]
            results_current = results_array[n_start:n_end]

            if any(results_current['e1'] == NO_RESULT_FLAG):
                logger.error('missing angles for %d %d' % (cid,gid))
                continue

            g1_true = truth_current['g1'][0]
            g2_true = truth_current['g2'][0]
            g1_mean = numpy.mean(results_current['e1'])
            g2_mean = numpy.mean(results_current['e2'])
            g1_stdm = numpy.std(results_current['e1'],ddof=1)
            g2_stdm = numpy.std(results_current['e2'],ddof=1)            

            
            shears_true_g1.append( g1_true )
            shears_true_g2.append( g2_true )
            shears_mean_g1.append( g1_mean )
            shears_mean_g2.append( g2_mean )
            shears_stdm_g1.append( g1_stdm )
            shears_stdm_g2.append( g2_stdm )

        n_valid_shears = len(shears_stdm_g2)
        logger.info('%d galaxy %d valid shears %d' % (i,cid,n_valid_shears))
        if n_valid_shears < 2:
            continue
            logger.error('not enough shears for galaxy %d %d : %d' % (i,cid,n_valid_shears))
        else:
            bias_g1 = numpy.array(shears_mean_g1,dtype=numpy.float64) - numpy.array(shears_true_g1,dtype=numpy.float64)
            bias_g2 = numpy.array(shears_mean_g2,dtype=numpy.float64) - numpy.array(shears_true_g2,dtype=numpy.float64)
            g1_true = numpy.array(shears_true_g1,dtype=numpy.float64) 
            g2_true = numpy.array(shears_true_g2,dtype=numpy.float64)
            bias_g1_stdm = numpy.ones(bias_g1.shape) 
            bias_g2_stdm = numpy.ones(bias_g2.shape)

            import pdb; pdb.set_trace()
            c1,m1,cov1 = _getLineFit(g1_true,bias_g1,bias_g1_stdm)
            c2,m2,cov2 = _getLineFit(g2_true,bias_g2,bias_g2_stdm)








    

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
    file_pickle = open(filename_pickle,'w')
    logger.info('saving')
    pickle.dump(results_array,file_pickle)
    file_pickle.close()
    # check if OK
    results_array = pickle.load(open(filename_pickle))
    logger.info('saved %s correctly, got %d rows' % (filename_pickle,results_array.shape[0]))
    file_pickle.close()


if __name__ == "__main__":

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
    getBiasForEachGal()
