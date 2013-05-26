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
import tabletools
import cPickle as pickle
from tablespec import *

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

def getBiasForEachGal():

    # get the results and truth
    results_array = tabletools.loadTable('results_array',args.filepath_results,dtype_table_results)
    truth_array   = tabletools.loadTable('truth_array',args.filepath_truth,dtype_table_truth)

    cosmos_ids = set(truth_array['id_cosmos'])
    n_cosmos_ids = len(cosmos_ids)
    n_per_cosmos = sum(truth_array['id_cosmos'] == truth_array[0]['id_cosmos'])
    print "n_per_cosmos" , n_per_cosmos
    n_shears = len(set(truth_array['id_shear']))
    n_angles = len(set(truth_array['id_angle']))

    logger.info('using set of %d cosmos galaxies ' % len(cosmos_ids) )

    results_stats = numpy.zeros(1,dtype=dtype_table_stats)

    n_gals_lost = 0

    for cid in range(n_cosmos_ids):

        shears_true_g1 = []
        shears_true_g2 = []
        shears_mean_g1 = []
        shears_mean_g2 = []
        shears_stdm_g1 = []
        shears_stdm_g2 = []
        rgp_list       = []
        snr_list       = []
        hlr_list       = []

        n_gal_start  = cid*n_per_cosmos
        zphot =  truth_array[n_gal_start]['zphot']
        id_cosmos = truth_array[n_gal_start]['id_cosmos']

        if cid % 100 == 0 : logger.info('passing %10d %10d %10d' % (cid,truth_array[n_gal_start]['id_unique'],results_array[n_gal_start]['identifier']))

        for gid in range(n_shears):

            n_start = cid*n_per_cosmos + gid     * n_angles
            n_end   = cid*n_per_cosmos + (gid+1) * n_angles

            truth_current   = truth_array[n_start:n_end]
            results_current = results_array[n_start:n_end]

            if any(results_current['e1'] == NO_RESULT_FLAG):
                logger.error('% 8d %10d -- missing angles for shear %d' % (cid,id_cosmos,gid))
                continue

            mean_g1 = numpy.mean(results_current['e1'])      
            mean_g2 = numpy.mean(results_current['e2'])      
            
            if numpy.isnan(mean_g1) or numpy.isnan(mean_g2):
                logger.error('% 8d %10d -- shear is nan for %d' % (cid,id_cosmos,gid))
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

        if n_valid_shears != n_shears:
            n_gals_lost+=1
            logger.error('% 8d %10d -- not enough shears for galaxy : %d , so far lost %d' % (cid,id_cosmos,n_valid_shears,n_gals_lost))
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

            results_stats_row = numpy.array([(cid,id_cosmos,zphot,m1,m2,m1_std,m2_std,c1,c2,c1_std,c2_std,hlr,rgp,snr)],dtype=dtype_table_stats)
            results_stats = numpy.concatenate((results_stats,results_stats_row))


        logger.debug('%8d galaxy %8d valid shears %d, so far got %8d valid gals' % (cid,id_cosmos,n_valid_shears,results_stats.shape[0]))

    # remove the first zero row
    results_stats=results_stats[1:]
    n_stats = len(results_stats)
    logger.info('results_stats has %d rows' % n_stats)

    file_pickle = open(args.filepath_stats,'w')
    pickle.dump(results_stats,file_pickle,protocol=2)
    file_pickle.close()
    logger.info('saved %s' % args.filepath_stats) 
    
def mergeResults():

    # truth_array   = loadTruthArray()
    truth_array = tabletools.loadTable('truth_array',args.filepath_truth,dtype_table_truth)
    n_gals_total = len(truth_array)
    n_gals_per_file = 640

    # get the wildcard for the files - join all files in a big array
    filecard_resutls = os.path.join(config['args'].dirpath_results,'results.nmb_main.real.*.cat')
    files = [os.path.join(config['args'].dirpath_results,'results.nmb_main.real.%05d.cat' % n) for n in range(0,n_gals_total,n_gals_per_file)]
    # files = glob.glob(filecard_resutls)
    # files.sort()
    
    n_files = len(files)
    n_meas = len(dtype_table_results['names'])
    logger.info('got %d file names' % n_files)

    # initialise the empty array
    global results_array
    results_array = numpy.zeros(n_gals_total,dtype=dtype_table_results)
    results_array['e1'] = NO_RESULT_FLAG 
    results_array['e2'] = NO_RESULT_FLAG
    results_array['time_taken'] = NO_RESULT_FLAG

    for fi,file_results in enumerate(files):

            if os.path.isfile(file_results):

                index_start = fi*n_gals_per_file
                index_end   = (fi+1)*n_gals_per_file
                results = numpy.loadtxt(file_results,dtype=dtype_table_results)
                n_gals_in_file = len(results)
                logger.info('%d file %s n_gals %d' % (fi,file_results,n_gals_in_file))
                if n_gals_in_file != n_gals_per_file:
                    logger.error('skipping file not enough galaxies')
                    continue
                else:

                    results_array[range(index_start,index_end)] = results
            else:
                logger.error('skipping file doesnt exist')
                continue

    logger.info('results  n %10d first %d last %d' % (len(results_array),results_array[0]['identifier'] ,results_array[-1]['identifier']))
    logger.info('truth    n %10d first %d last %d' % (len(truth_array),truth_array[0]['id_unique']      ,truth_array[-1]['id_unique']))

    n_matches = 0
    for ri in range(len(results_array)):
        if results_array[ri]['identifier'] == truth_array[ri]['id_unique']: n_matches+=1
    logger.info('number of matches %s' % n_matches)
    
    results_array['identifier'] = truth_array['id_unique']
    
    n_matches = 0
    for ri in range(len(results_array)):
        if results_array[ri]['identifier'] == truth_array[ri]['id_unique']: n_matches+=1
    logger.info('number of matches %s' % n_matches)

    # save file pickle - use global filename_pickle
    tabletools.saveTable(args.filepath_results,results_array)
    logger.info('results saved %s correctly, got %d rows' % (args.filepath_results,len(results_array)))
    
    tabletools.saveTable(args.filepath_truth,truth_array)
    logger.info('truth saved %s correctly, got %d rows' % (args.filepath_truth,len(truth_array)))

    logger.info('truth   n %10d first %d last %d' % (len(results_array),results_array[0]['identifier'],results_array[-1]['identifier']))
    logger.info('results n %10d first %d last %d' % (len(truth_array),truth_array[0]['id_unique'],truth_array[-1]['id_unique']))



def createBFITsample():

    # load the truth
    truth_array = tabletools.loadTable('truth_array',args.filepath_truth,dtype=dtype_table_truth)
    logger.info('truth array have %d galaxies' % len(truth_array))

    
    # load the results
    results_array = tabletools.loadTable('results_array',args.filepath_results,dtype=dtype_table_results)
    logger.info('results array have %d galaxies' % len(results_array))

    # load the stats - they contain the galaxies which passed
    stats_array = tabletools.loadTable('stats_sarray',args.filepath_stats,dtype=dtype_table_stats)
    logger.info('stats array have %d galaxies' % len(stats_array))

    results_array_bfit = numpy.zeros(1,dtype=dtype_table_results)
    truth_array_bfit = numpy.zeros(1,dtype=dtype_table_truth)

    for ig,g in enumerate(stats_array['cosmos_id']):
        select = truth_array['id_cosmos'] == g

        results_rows = results_array[select] 
        truth_rows = truth_array[select]

        if ig % 100 == 0 : 
            logger.info('passing galaxy %10d results id %10d truth id %10d' % (ig,results_rows['identifier'][0],truth_rows['id_unique'][0]))
            
        results_array_bfit = numpy.append(results_array_bfit,results_rows)
        truth_array_bfit = numpy.append(truth_array_bfit,truth_rows)

    # remove the last one
    results_array_bfit = results_array_bfit[1:]
    truth_array_bfit = truth_array_bfit[1:]
    n_gals = len(truth_array_bfit)

    logger.info('number of galaxies in bfit sample %d' % n_gals)

    n_gals = len(stats_array)

    filename_results_bfit = args.filepath_results.replace('results','bfit').replace('pp','fits')
    filename_truth_bfit = 'truth.%d.fits' % n_gals

    fits_results = tabletools.getFITSTable(results_array_bfit)
    fits_truth   = tabletools.getFITSTable(truth_array_bfit)
    n_gals_fits = fits_truth[1].data.shape
    logger.info('got fits table from numpy, with %d rows' % n_gals_fits)

    fits_results.writeto(filename_results_bfit,clobber=True)
    fits_truth.writeto(filename_truth_bfit,clobber=True)
    logger.info('saved %s %s' % (filename_truth_bfit,filename_results_bfit))

    # check

    results_array_bfit_loaded = pyfits.open(filename_results_bfit)
    n_loaded = len(results_array_bfit_loaded[1].data)
    logger.info('loaded n_gals %d' % n_loaded)
    logger.info('first ids %10d %10d' % (results_array_bfit_loaded[1].data[ 0]['identifier'] , results_array_bfit['identifier'][0]))
    logger.info('last  ids %10d %10d' % (results_array_bfit_loaded[1].data[-1]['identifier'] , results_array_bfit['identifier'][-1]))
    logger.info('first e1 % f % f' % (results_array_bfit_loaded[1].data['e1'][0]  , results_array_bfit['e1'][0]  ))
    logger.info('last  e1 % f % f' % (results_array_bfit_loaded[1].data['e1'][-1] , results_array_bfit['e1'][-1] ))




def main():

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('command', type=str, help='what to do?')
    parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    parser.add_argument('--filepath_truth', type=str, default='truth.26000.pp', help='truth file for the run, overrides the config file (by default is taken from yaml file)')
    parser.add_argument('--filepath_stats', type=str, default='stats.nmb_main.real.pp', help='stats file')
    parser.add_argument('--filepath_results', type=str, default='results.nmb_main.real.pp', help='results file')
    parser.add_argument('--dirpath_results', type=str, default='results', help='where the results files are')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    global logger , config , args

    args = parser.parse_args()
    args.name_config = os.path.basename(args.filepath_config).replace('.yaml','')

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("nmb_main") 
    logger.setLevel(logging_level)
    
    # load the configuration file
    config = yaml.load(open(args.filepath_config,'r')) 
    
    eval(args.command + '()')
    # mergeResutls()
    # getBiasForEachGal()
    # createBFITsample()

if __name__ == "__main__":

    main()