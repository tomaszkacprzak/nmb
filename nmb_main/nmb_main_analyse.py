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

def loadTable(table_name,filepath,dtype):

    try:
        table = eval(table_name)
    except:

        logger.info('loading %s' % filepath)
        if filepath.split('.')[-1] == 'pp':
                file_pickle = open(filepath)
                table = pickle.load(file_pickle)
                file_pickle.close()
        else:
                table = numpy.loadtxt(filepath,dtype=dtype)
        
    else:
        logger.info('using preloaded array %s' % table_name)
    
    logger.info('loaded %s correctly, got %d rows' % (filepath,len(table)))

    globals()[table_name] = table
    return table



def loadResultsArray():
    """
    Loads the resutls array from file global filename_results_pickle, unless there is a lobal results_array.
    @return results_array
    """

    global results_array
    
    try:
        results_array
    except:
        logger.info('loading %s' % args.filepath_results)
        if args.filepath_results.split('.')[-1] == 'pp':
                file_pickle = open(args.filepath_results)
                results_array = pickle.load(file_pickle)
                file_pickle.close()
        else:
                results_array = numpy.loadtxt(args.filepath_results,dtype=dtype_table_results)
        
    else:
        logger.info('using preloaded results_array')
    
    logger.info('loaded %s correctly, got %d rows' % (args.filepath_results,results_array.shape[0]))

    return results_array

def loadTruthArray():
    """
    Loads the truth array from file global filename_truth_pickle, unless there is a global truth_array.
    @return truth_array
    """

    global truth_array
    
    try:
        truth_array
    except:
        logger.info('loading %s' % args.filepath_truth)
        if args.filepath_truth.split('.')[-1] == 'pp':
            file_pickle = open(args.filepath_truth)
            truth_array = pickle.load(file_pickle)
            file_pickle.close()
        else:
            truth_array = numpy.loadtxt(args.filepath_truth,dtype=dtype_table_truth)

    else:
        logger.info('using preloaded truth_array')
    
    logger.info('loaded %s correctly, got %d rows' % (args.filepath_truth,truth_array.shape[0]))

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

        n_gal_start  = i*n_per_cosmos

        if i % 100 == 0 : logger.info('passing %d %d %d %d' % (i,cid,truth_array[n_gal_start]['id_unique'],results_array[n_gal_start]['identifier']))

        for gid in range(n_shears):

            n_start = i*n_per_cosmos + gid     * n_angles
            n_end   = i*n_per_cosmos + (gid+1) * n_angles

            zphot =  truth_array[n_start]['zphot']

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
        if n_valid_shears != n_shears:
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

            # dtype_table_stats =  { 'names'   : ['index', 'cosmos_id','zphot','m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std' , 'hlr' , 'rgp' , 'snr'],
            results_stats_row = numpy.array([(i,cid,zphot,m1,m2,m1_std,m2_std,c1,c2,c1_std,c2_std,hlr,rgp,snr)],dtype=dtype_table_stats)
            results_stats = numpy.concatenate((results_stats,results_stats_row))


        logger.debug('%8d galaxy %8d valid shears %d, so far got %8d valid gals' % (i,cid,n_valid_shears,results_stats.shape[0]))

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
    truth_array = loadTable('truth_array',args.filepath_truth,dtype_table_truth)

    # get the wildcard for the files - join all files in a big array
    filecard_resutls = os.path.join(config['args'].dirpath_results,'results.nmb_main.real.*.cat')
    files = [os.path.join(config['args'].dirpath_results,'results.nmb_main.real.%05d.cat' % n) for n in range(0,n_gals_total,n_gals_per_file)]
    # files = glob.glob(filecard_resutls)
    # files.sort()
    
    n_gals_total = len(truth_array)
    n_gals_per_file = 640
    n_files = len(files)
    n_meas = len(dtype_table_results['names'])
    logger.info('got %d file names' % n_files)

    # initialise the empty array
    global results_array
    results_array = numpy.zeros(n_gals_total,dtype=dtype_table_results)
    results_array['identifier'] = truth_array['id_unique']
    results_array['e1'] = NO_RESULT_FLAG 
    results_array['e2'] = NO_RESULT_FLAG
    results_array['time_taken'] = NO_RESULT_FLAG

    for fi,file_results in enumerate(files):

            index_start = fi*n_gals_per_file
            index_end   = (fi+1)*n_gals_per_file
            results = numpy.loadtxt(file_results,dtype=dtype_table_results)
            n_gals_in_file = len(results)
            logger.info('%d file %s n_gals %d' % (fi,file_results,n_gals_in_file))
            if n_gals_in_file != n_gals_per_file:
                continue
            else:
                results_array[index_start,index_end] = results











    # results_all = numpy.zeros(1,dtype=dtype_table_results)

    # filename_results_all = 'results_all.pp'
    # if not os.path.isfile(filename_results_all):

    #     # loop over files
    #     for fi,file_results in enumerate(files):

    #         results = numpy.loadtxt(file_results,dtype=dtype_table_results)
    #         n_gals_in_file = len(results)
    #         results_all = numpy.append(results_all,results)      
            
    #         if fi % 100 == 0 : logger.info('%4d loaded file %50s with %4d lines, results got so far %d' %(fi,file_results,n_gals_in_file,len(results_all)))

    #     pickle.dump(results_all,open('results_all.pp','w'))
    # else:
    #     results_all = pickle.load(open(filename_results_all,'r'))

    # # remove the ones that are in the results and not in the truth
    # sr = set(results_all['identifier'])
    # st = set(truth_array['id_unique'])
    # diff_set = sr.difference(st)
    # logger.info('found %d elements in diff set' % len(diff_set))

    # logger.info('results_all has %d elements' % len(results_all))
    # for di,ds in enumerate(diff_set):

    #     select = numpy.nonzero(results_all['identifier'] == ds)     
    #     results_all = numpy.delete(results_all,select)
    #     logger.info('%d diff set id %d removing %d elements' % (di,ds,len(select)))

    # logger.info('removed diff set - results_all has %d elements' % len(results_all))

    # # remove first empty row
    # truth_array_sorted = truth_array[numpy.argsort(truth_array['id_unique'])]
    # results_all_sorted = results_all[numpy.argsort(results_all['identifier'])]
    # logger.info('results  n %10d first %d last %d' % (len(results_array),results_all_sorted[0]['identifier'],results_all_sorted[-1]['identifier']))
    # logger.info('truth    n %10d first %d last %d' % (len(truth_array),truth_array_sorted[0]['id_unique'],truth_array_sorted[-1]['id_unique']))

    # n_gals_avail = len(results_all_sorted)

    # n_range = 2000;
    # for i,idu in enumerate(truth_array_sorted['id_unique']):

    #     n_start = max(0,i-n_range)
    #     n_end = min(n_gals_avail,i+n_range)
    #     select = numpy.nonzero(results_all_sorted[range(n_start,n_end)]['identifier'] == idu)
    #     n_found = sum(select)
    #     if n_found == 1: 
    #         row = results_all_sorted[select]
    #         results_array[i] = row
    #         print row['identifier'] , len(row)
    #     elif n_found == 0:
    #         continue
    #     else:
    #         rows_found = results_all_sorted[select]
    #         raise ValueError('n_found %d value %d' % (n_found,rows_found[0]['identifier']))

    #     if i % 100 == 0:  logger.info('finished %d, found %d' % (i,len(select)))


    # logger.info('found %d results' % len(results_all))      

    # j=0
    # # match the catalogs
    # for i,idu in enumerate(truth_array_sorted['id_unique']):

    #     if i % 10000 == 0 : logger.info('passing %d %d' % (i,j))

    #     if idu == results_all_sorted[j]['identifier']:
    #         results_array[i] = results_all_sorted[j]
    #         j+=1
    #     else:
    #         logger.debug('%5d %5d %10d %10d not found' % (i,j,idu,results_all_sorted[j]['identifier']))
    #         continue

    logger.info('finished %d %d' % (i,j))

    # save file pickle - use global filename_pickle
    saveTable(args.filepath_results,results_array)
    logger.info('results saved %s correctly, got %d rows' % (args.filepath_results,len(results_array.shape)))
    
    saveTable(args.filepath_truth,truth_array)
    logger.info('truth saved %s correctly, got %d rows' % (filename_truth_pickle,len(truth_array.shape)))

    logger.info('truth   n %10d first %d last %d' % (len(results_array),results_array[0]['identifier'],results_array[-1]['identifier']))
    logger.info('results n %10d first %d last %d' % (len(truth_array),truth_array[0]['id_unique'],truth_array[-1]['id_unique']))

def saveTable(filepath,table):

    if filepath.split('.')[-1] == 'pp':
        file_pickle = open(filepath,'w')
        pickle.dump(table,filepath,protocol=2)
        file_pickle.close()
    else:
        header = '# ' + ' '.join(table.dtype.names)
        numpy.savetxt(filepath,table,header=header)

    logger.info('truth saved %s correctly, got %d rows' % (filepath,len(table.shape)))



def createBFITsample():

    # load the stats - they contain the galaxies which passed
    file_stats_pickle = open(filename_stats_pickle)
    results_stats = pickle.load(file_stats_pickle)
    logger.info('results stats have %d galaxies' % len(results_stats))

    # load the resutls
    results_array = loadResultsArray()
    logger.info('results array have %d galaxies' % len(results_array))

    truth_array = loadTruthArray()
    logger.info('truth array have %d galaxies' % len(truth_array))

    results_array_bfit = numpy.zeros(1,dtype=dtype_table_results)

    for ig,g in enumerate(results_array):
        if ig % 100 == 0 : logger.info('passing galaxy %d' % ig)
        if truth_array[ig]['id_cosmos'] in results_stats['cosmos_id']:
                import pdb; pdb.set_trace()
                results_row = numpy.array(results_array[ig],dtype=dtype_table_results)

    # remove the last one
    results_array_bfit = results_array_bfit[1:]
    n_gals = len(results_array_bfit)

    logger.info('number of galaxies in bfit sample %d' % n_gals)

    fits = getFITSTable(results_array_bfit)
    n_gals_fits = fits[1].data.shape
    logger.info('got fits table from numpy, with %d rows' % n_gals_fits)

def getFITSTable(numpy_array):

    formats = { numpy.dtype('int32') : 'J' , numpy.dtype('float32') : 'D' }

    cols = []

    for i,col_name in enumerate(numpy_array.dtype.names):
        col_type = numpy_array.dtype[i]
        col_fmt = formats[col_type]
        print col_fmt
        col = pyfits.Column(name=col_name,format=col_fmt,array=numpy_array[col_name])
        cols.append(col)

    print cols
    tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
    import pdb; pdb.set_trace()

    hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([hdu, tbhdu])



def main():

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('command', type=str, help='what to do?')
    parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    parser.add_argument('--filepath_truth', type=str, default='truth.26000.cat', help='truth file for the run, overrides the config file (by default is taken from yaml file)')
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
    # store the args in config so it's easier to use them
    config['args'] = args    

    eval(args.command + '()')
    # mergeResutls()
    # getBiasForEachGal()
    # createBFITsample()

if __name__ == "__main__":

    main()