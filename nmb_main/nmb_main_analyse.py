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
import glob
import tabletools
import cPickle as pickle
from tablespec import *

default_logger = logging.getLogger('nmb_main_analyse')
default_logger.setLevel(logging.WARNING)

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

def getTotalBias():

    dtype_table_results_use = dtype_table_results2  

    # get the results and truth
    results_array = tabletools.loadTable(table_name='results_array',filepath=args.filepath_results,dtype=dtype_table_results_use)
    truth_array   = tabletools.loadTable(table_name='truth_array'  ,filepath=args.filepath_truth  ,dtype=dtype_table_truth)    
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)   

    mc_results_bfit  = getBiasForResults(results_real_noisy,truth_array,logger)
    mc_results_real = getBiasForResults(results_bfit_noisy,truth_array,logger)
    print 'real %2.6f %2.6f +/- %2.6f %2.6f' % (mc_results_real['m1'],mc_results_real['m2'],mc_results_real['m1_std'],mc_results_real['m2_std'])    
    print 'bfit %2.6f %2.6f +/- %2.6f %2.6f' % (mc_results_bfit['m1'],mc_results_bfit['m2'],mc_results_bfit['m1_std'],mc_results_bfit['m2_std'])    
    




def getBiasForResults(results_array,truth_array,n_means=5,bin_param='test',bin_id=0,logger=default_logger):

    

    # see how many shears and angles we have
    n_shears = len(set(truth_array['id_shear']))
    n_angles = len(set(truth_array['id_angle']))

    true_g1_list = [truth_array['g1'][ truth_array['id_shear'][0:n_shears*n_angles] == idg ][0]  for idg in range(n_shears)]
    true_g2_list = [truth_array['g2'][ truth_array['id_shear'][0:n_shears*n_angles] == idg ][0]  for idg in range(n_shears)]


    # those will contain len(results_array)/n_gals_per_mean points
    shears_true_g1 = []  
    shears_true_g2 = []  
    shears_mean_g1 = []  
    shears_mean_g2 = []  
    shears_stdm_g1 = []  
    shears_stdm_g2 = []  
    shears_stdv_g1 = []  
    shears_stdv_g2 = []  
    shears_n_gals  = [] 

    # loop over shears
    for sid in range(n_shears):

        # select shear
        if 'id_cosmos' in results_array.dtype.names:
            select = getShearIDfromUniqueID(results_array['id_unique'])  == sid
        elif 'identifier' in  results_array.dtype.names:
            select = getShearIDfromUniqueID(results_array['identifier'])  == sid

        results_sid = results_array[select]
        n_gals = len(results_sid)
        n_gals_per_mean = n_gals/n_means
        # n_means = numpy.ceil(float(n_gals) / float(n_gals_per_mean)).astype(numpy.int64) 
        logger.debug('shear %d n_gals %d n_means %d' % (sid,n_gals,n_means))



        for idm in range(n_means):          

            i_start = idm    *n_gals_per_mean
            i_end   = (idm+1)*n_gals_per_mean if (idm+1)*n_gals_per_mean < n_gals else n_gals


            results_sid_part = results_sid[i_start:i_end]         
            # remove all the fields with errors
            select = numpy.logical_and(results_sid_part['e1'] < 1.,results_sid_part['e2'] < 1.)
            results_sid_use = results_sid_part[select]
            n_gals_use = len(results_sid_use)

            mean_g1 = numpy.mean(results_sid_use['e1'])
            mean_g2 = numpy.mean(results_sid_use['e2'])
            stdv_g1 = numpy.std(results_sid_use['e1'],ddof=1) 
            stdv_g2 = numpy.std(results_sid_use['e2'],ddof=1) 
            stdm_g1 = stdv_g1 / numpy.sqrt(n_gals_use)
            stdm_g2 = stdv_g2 / numpy.sqrt(n_gals_use)

            logger.debug('%10d %10d \t% 2.4f\t% 2.4f\t% 2.4f\t% 2.4f\t% 2.4f\t% 2.4f \t using %10d / %10d, min/max e1/e2 = %2.2f %2.2f %2.2f %2.2f' % (i_start,i_end,mean_g1,mean_g2,stdm_g1,stdm_g2,stdv_g1,stdv_g2,len(results_sid_use),len(results_sid_part),results_sid_use['e1'].min(),results_sid_use['e1'].max(),results_sid_use['e2'].min(),results_sid_use['e2'].max()))         

            shears_true_g1.append( true_g1_list[sid] )
            shears_true_g2.append( true_g2_list[sid] )
            shears_mean_g1.append( mean_g1 )
            shears_mean_g2.append( mean_g2 )
            shears_stdm_g1.append( stdm_g1 )
            shears_stdm_g2.append( stdm_g2 )
            shears_stdv_g1.append( stdv_g1 )
            shears_stdv_g2.append( stdv_g2 )
            shears_n_gals.append(n_gals_use)

    tg1 = numpy.array(shears_true_g1)
    mg1 = numpy.array(shears_mean_g1)
    sg1 = numpy.array(shears_stdm_g1)
    bg1 = mg1 - tg1;
    c1,m1,cov1 = _getLineFit(tg1,bg1,sg1)
    lg1 = tg1*m1 + c1
    zg1 = numpy.ones(sg1.shape)*numpy.std(lg1 - bg1,ddof=1)
    c1,m1,cov1 = _getLineFit(tg1,bg1,zg1)
    std_m1 = numpy.sqrt(cov1[1,1])
    std_c1 = numpy.sqrt(cov1[0,0])

    tg2 = numpy.array(shears_true_g2)
    mg2 = numpy.array(shears_mean_g2)
    sg2 = numpy.array(shears_stdm_g2)
    bg2 = mg2 - tg2;
    c2,m2,cov2 = _getLineFit(tg2,bg2,sg2)
    lg2 = tg2*m2 + c2
    zg2 = numpy.ones(sg2.shape)*numpy.std(lg2 - bg2,ddof=1)
    c2,m2,cov2 = _getLineFit(tg2,bg2,zg2)
    std_m2 = numpy.sqrt(cov2[1,1])
    std_c2 = numpy.sqrt(cov2[0,0])    

    if logger!=None:
        if logger.level <= logging.DEBUG:
            import pylab
            
            pylab.figure()
            pylab.errorbar(tg1,bg1,yerr=sg1,fmt='r.')
            pylab.errorbar(tg1,bg1,yerr=zg1,fmt='m.')
            indices = tg1.argsort()
            pylab.plot(tg1[indices],lg1[indices],'r-')
            # pylab.plot(tg1[indices],tg1[indices],'-k')
            
            pylab.errorbar(tg2,bg2,yerr=sg2,fmt='b.')
            pylab.errorbar(tg2,bg2,yerr=zg2,fmt='c.')
            indices = tg2.argsort()
            pylab.plot(tg2[indices],lg2[indices],'b-')
            # pylab.plot(tg1[indices],tg1[indices],'-k')

            filename_fig = 'debug/fig.bins.%s.%02d.png' % (bin_param,bin_id)

            pylab.xlabel('g_true')
            pylab.ylabel('g_est - g_true')
            pylab.title('bin %d in %s' % (bin_id,bin_param))

            pylab.savefig(filename_fig)
            pylab.close()
            logger.info('saved file %s' % filename_fig)


    
    logger.info('m1 = % 2.4f \t +/- % 2.4f' % ( m1, std_m1))
    logger.info('m2 = % 2.4f \t +/- % 2.4f' % ( m2, std_m2))
    logger.info('c1 = % 2.4f \t +/- % 2.4f' % ( c1, std_c1))
    logger.info('c2 = % 2.4f \t +/- % 2.4f' % ( c2, std_c2))

    # dtype_table_binstats = { 'names'   : ['bin_id', 'bin_value' , 'bin_ngals' ,'m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std'],
    mc_results = numpy.array([(0,0,len(results_array),m1,m2,std_m1,std_m2,c1,c2,std_c1,std_c2)],dtype=dtype_table_binstats)

    return mc_results


def getShearIDfromUniqueID(ids_unique):

    # enforce numpy array
    ids = numpy.array(ids_unique)
    return ((ids % 10000) // 100).astype(numpy.int64)

def getCosmosIDfromUniqeID(ids_unique):
    ids = numpy.array(ids_unique)
    return ((ids / 10000)).astype(numpy.int64)


def getBiasForEachGal():

    # get the results and truth
    results_array = tabletools.loadTable(table_name='results_array',filepath=args.filepath_results,dtype=dtype_table_results)
    truth_array   = tabletools.loadTable(table_name='truth_array'  ,filepath=args.filepath_truth  ,dtype=dtype_table_truth)

    cosmos_ids = set(truth_array['id_cosmos'])
    n_cosmos_ids = len(cosmos_ids)
    n_per_cosmos = sum(truth_array['id_cosmos'] == truth_array[0]['id_cosmos'])
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
    """
    @brief read all the results files, and then merge them together and save it in a fits binary table.
    """

    n_reps = args.n_reps 

    # truth_array   = loadTruthArray()
    truth_array = tabletools.loadTable(table_name='truth_array',filepath=args.filepath_truth,dtype=dtype_table_truth)
    # n_gals_total = len(truth_array)
    n_gals_total = config['settings']['n_images']*n_reps
    logger.info('using %d repetitions of the all config galaxies, total %d gals' % (n_reps,n_gals_total) )
    n_gals_per_file = args.n_gals_per_file

    # get the wildcard for the files - join all files in a big array
    filecard = 'results.%s' % config['name']  + '.%012d.cat'
    files1 = [os.path.join(args.dirpath_results,filecard % n) for n in range(0,n_gals_total,n_gals_per_file)]
    filecard = 'results.%s' % config['name']  + '.%05d.cat'
    files2 = [os.path.join(args.dirpath_results,filecard % n) for n in range(0,n_gals_total,n_gals_per_file)]

    # do a small test to get correct format

    if os.path.isfile(files1[0]) : file_results = files1[0]
    elif os.path.isfile(files2[0]) : file_results = files2[0]
    results = numpy.loadtxt(file_results)
    if results.shape[1] == 23:
        dtype_table_results_use = dtype_table_results2
        id_unique_field_results = 'id_unique'
        logger.info('using dtype_table_results2 with %d columns' % results.shape[1])
    elif results.shape[1] == 20:
        dtype_table_results_use = dtype_table_results
        id_unique_field_results = 'identifier'
        logger.info('using dtype_table_results with %d columns' % results.shape[1])
    else:
        logger.error('results have unknown column structure')

    # get the dimensionalty of the datas - doesn't matter if we use files1 or files2 for that
    n_files = len(files1)
    n_meas = len(dtype_table_results['names'])
    logger.info('got %d file names' % n_files)

    # initialise the empty array
    global results_array
    results_array = numpy.zeros(n_gals_total,dtype=dtype_table_results_use)
    results_array['e1'] = NO_RESULT_FLAG 
    results_array['e2'] = NO_RESULT_FLAG
    results_array['time_taken'] = NO_RESULT_FLAG

    n_files_found = 0
    n_files_valid = 0

    for fi in range(len(files1)):

            if os.path.isfile(files1[fi]) or os.path.isfile(files2[fi]):

                n_files_found += 1
                
                if os.path.isfile(files1[fi]) : file_results = files1[fi]
                elif os.path.isfile(files2[fi]) : file_results = files2[fi]
    
                # locate the part of the total results_array that this file corresponds to
                index_start = fi*n_gals_per_file
                index_end   = (fi+1)*n_gals_per_file
                results = numpy.loadtxt(file_results,dtype=dtype_table_results_use)
                n_gals_in_file = len(results)
                if n_gals_in_file != n_gals_per_file:
                    logger.error('%5d file %s not enough galaxies ----- skipping, so far found %d files, missing %d, accepted %d , invalid %d' % (fi,file_results,n_files_found,fi-n_files_found+1,n_files_valid,fi-n_files_valid+1))
                else:
                    n_files_valid += 1
                    logger.info('%5d file %s n_gals %d' % (fi,file_results,n_gals_in_file))
                    results_array[range(index_start,index_end)] = results
            else:
                file_results = files1[fi]
                logger.error('%5d file %s file not found      ----- skipping, so far found %d files, missing %d, accepted %d , invalid %d' % (fi,file_results,n_files_found,fi-n_files_found+1,n_files_valid,fi-n_files_valid+1))
                

    logger.info('results  n %10d first %d last %d' % (len(results_array),results_array[0][id_unique_field_results] ,results_array[-1][id_unique_field_results]))
    logger.info('truth    n %10d first %d last %d' % (len(truth_array),truth_array[0][id_unique_field_results]     ,truth_array[-1][id_unique_field_results]))
    logger.info('getting number of matches ...')
    
    # get the number of matches before the assignment of IDS
    n_matches = 0
    n_gals_in_truth = len(truth_array)
    for ri in range(len(results_array)):
        # roll the ri_truth to be always within the len of truth_array, as results may have more reps
        if ri >= n_gals_in_truth: 
            ri_truth = ri % n_gals_in_truth
        else: ri_truth = ri
        if results_array[ri][id_unique_field_results] == truth_array[ri_truth]['id_unique']: n_matches+=1
    logger.info('number of matches %s' % n_matches)
    
    # assign the ids
    # import pdb; pdb.set_trace()
    truth_id_extended_array = numpy.concatenate([truth_array['id_unique']]*n_reps)
    results_array[id_unique_field_results] = truth_id_extended_array
    
    # get the number of matches after the assignment of IDS  
    n_matches = 0
    for ri in range(len(results_array)):
        # roll the ri_truth to be always within the len of truth_array, as results may have more reps
        if ri >= len(truth_array): 
            ri_truth = ri % len(truth_array)
        else: ri_truth = ri
        if results_array[ri][id_unique_field_results] == truth_array[ri_truth]['id_unique']: n_matches+=1    
    logger.info('number of matches %s' % n_matches)
    logger.info('saving tables ...')
    tabletools.saveTable(args.filepath_results,results_array)
    logger.info('results saved %s correctly, got %d rows' % (args.filepath_results,len(results_array)))
    
    # tabletools.saveTable(args.filepath_truth,truth_array)
    # logger.info('truth saved %s correctly, got %d rows' % (args.filepath_truth,len(truth_array)))

    logger.info('results   n %10d first %d last %d' % (len(results_array),results_array[0][id_unique_field_results],results_array[-1][id_unique_field_results]))
    logger.info('truth     n %10d first %d last %d' % (len(truth_array),truth_array[0]['id_unique'],truth_array[-1]['id_unique']))

def createBFITsample():

    # load the truth
    truth_array = tabletools.loadTable(table_name='truth_array',filepath=args.filepath_truth,dtype=dtype_table_truth)
    logger.info('truth array have %d galaxies' % len(truth_array))

    
    # load the results
    results_array = tabletools.loadTable(table_name='results_array',filepath=args.filepath_results,dtype=dtype_table_results)
    logger.info('results array have %d galaxies' % len(results_array))

    # load the stats - they contain the galaxies which passed
    stats_array = tabletools.loadTable(table_name='stats_sarray',filepath=args.filepath_stats,dtype=dtype_table_stats)
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

def selectByIDs(ids,results_array,truth_array,logger=default_logger):
    """
    @brief selects a subset of results which correspond to the ids
    @param ids a list of integer ids for cosmos galaxies
    @param results_array results array to be selected from
    @param truth_array results array to be selected from
    @return results_bin results,truth_bin,indices_bin
    """

    # get the counts for various things
    cosmos_ids = set(truth_array['id_cosmos'])
    n_cosmos_ids = len(cosmos_ids)
    n_per_cosmos = sum(truth_array['id_cosmos'] == truth_array[0]['id_cosmos'])
    n_shears = len(set(truth_array['id_shear']))
    n_angles = len(set(truth_array['id_angle']))

    # get lenth of results_array
    n_measurements = len(results_array)

    # get the cosmos ids, but unique in the same order
    select = range(0,n_measurements,n_per_cosmos)
    if 'id_unique' in results_array.dtype.names:    ids_cosmos_unique = results_array[select]['id_unique']       
    elif 'identifier' in results_array.dtype.names: ids_cosmos_unique = results_array[select]['identifier']

    ids_cosmos_unique = getCosmosIDfromUniqeID(ids_cosmos_unique)

    # initialise a list to hold the indices
    indices_bin = []

    # loop over all cosmos ids in this bin
    for j,current_id in enumerate(ids):

        # get the index of the first occurence of this cosmos id in the ids_cosmos_unique
        id_indices_in_truth = numpy.nonzero(ids_cosmos_unique == current_id)[0]

        # logger.debug('selecting by ID %10d, found %d places in results list: %s' % (current_id,len(id_indices_in_truth),str(id_indices_in_truth)))

        for id_index_in_truth in id_indices_in_truth:

            # if found - append indices_bins with indices of that cosmos id in the results
            
            id_index_in_results_start = id_index_in_truth    *n_per_cosmos
            id_index_in_results_end   = (id_index_in_truth+1)*n_per_cosmos
            indices_bin.extend(range(id_index_in_results_start,id_index_in_results_end))
            

    # add the results
    if len(results_array) == len(truth_array):
        results_bin = results_array[indices_bin]
        truth_bin = truth_array[indices_bin]
    else :
        n_reps = len(results_array)/len(truth_array)
        logger.debug('n_reps = %d' % n_reps)
        results_bin = results_array[indices_bin]
        truth_array_rep = numpy.concatenate([truth_array]*n_reps)
        truth_bin = truth_array_rep[indices_bin] 

    logger.debug('%5d unique IDs given, found %5d entries in results_array' % (len(ids),len(results_bin)))

    return results_bin,truth_bin,indices_bin


def main():

    description = 'Plots for the results of nmb_main.'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('command', type=str, help='what to do?')
    parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    parser.add_argument('--filepath_truth', type=str, help='truth file for the run, overrides the config file (by default is taken from yaml file)')
    parser.add_argument('--filepath_stats', type=str, default='auto', help='stats file')
    parser.add_argument('--filepath_results', type=str, default='auto', help='results file')
    parser.add_argument('--dirpath_results', type=str, default='results', help='where the results files are')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('--n_gals_per_file', type=int, default=640 , help='number of galaxies in one results file from legion')
    parser.add_argument('--n_reps', type=int, default=1 , help='number of times the config file was ran')
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
    config['name'] = os.path.basename(args.filepath_config).replace('.yaml','')

    # get the filenames of stats and results using config name
    if args.filepath_results == 'auto':
        args.filepath_results = 'results.%s.fits' % config['name'] 
    if args.filepath_stats == 'auto':
        args.filepath_stats = 'stats.%s.fits' % config['name'] 

    
    eval(args.command + '()')
    # mergeResults()
    # getBiasForEachGal()
    # createBFITsample()

if __name__ == "__main__":

    main()