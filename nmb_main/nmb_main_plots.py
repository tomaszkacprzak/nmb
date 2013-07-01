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
import tabletools
from tablespec import *
import nmb_main_analyse as analyse

filepath_stats              = 'stats.nmb_main.real.pp'
filepath_acs_join_stats     = 'acs_joins_stats.pp'
filepath_truth_25880        = 'truth.25880.fits'
filepath_truth_26000        = 'truth.26000.pp'
filepath_acs                = 'cosmos_acs_shera_may2011.fits.gz'
filepath_results_real       = 'results.nmb_main.real.pp' 
# filepath_results_real_noisy = 'results.nmb_main.real.noisy.fits'
filepath_results_real_noisy = 'results.nmb_main.real.noisy.rep3.fits'
# filepath_results_bfit_noisy = 'results.nmb_main.bfit.noisy.fits'
filepath_results_bfit_noisy = 'results.nmb_main.bfit.noisy.rep2.fits'

NO_RESULT_FLAG = 666


bins_redshift = [0 , 0.35, 0.6, 0.8, 1.1 , 1.5]    
bins_size     = [1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8, 1.9 , 2.0]
bins_size_nocut= [1.0, 1.1 , 1.2, 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8, 1.9 , 2.0]
bins_modd     = [0 , 5 ,  9 , 20, 32]
bins_snr = numpy.logspace(2,4,20)
bins_hlr = numpy.linspace(0,5,5)

def getColorMap(n_colors):
    import colorsys
    HSV_tuples = [(x*1.0/n_colors, 0.75, 0.75) for x in range(n_colors)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def getTableACSjoinStats():


    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)
    results_modd = []
    
    logger.info('getting modd')
    for i in range(len(stats_array)):
        select = numpy.nonzero(acs_array['IDENT']==stats_array['cosmos_id'][i])
        results_modd.append( acs_array[select] )

    results_modd = numpy.concatenate(results_modd)
    logger.info('got modd')
    tabletools.saveTable(filepath_acs_join_stats,results_modd)
    return results_modd

def plotBiasForBins():

    req1_m = 0.02
    req2_m = 0.004

# redshift ------------------------------------------------------------------------------------------------------------

    binstats_real_noisy = tabletools.loadTable(filepath='bins.redshift.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.redshift.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.redshift.real.clear.cat',dtype=dtype_table_binstats)

    weights = binstats_real_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_real_noisy['bin_ngals']))
    logger.info( 'm1 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m1'] , weights ))
    logger.info( 'm2 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m2'] , weights ))

    weights = binstats_bfit_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_bfit_noisy['bin_ngals']))
    logger.info( 'm1 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m1'] , weights ))
    logger.info( 'm2 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m2'] , weights ))


    pylab.figure()
    pylab.clf()

    bins_centered = _binCentersSameLen(binstats_real_clear['bin_value'])

    pylab.errorbar(bins_centered,binstats_real_clear['m1'],yerr=binstats_real_clear['m1_std'], fmt='-y+', label = 'm1 model bias')
    pylab.errorbar(bins_centered,binstats_real_clear['m2'],yerr=binstats_real_clear['m2_std'], fmt='-g', label = 'm2 model bias')


    pylab.errorbar(bins_centered,binstats_bfit_noisy['m1'],yerr=binstats_bfit_noisy['m1_std'], fmt='-r+', label = 'm1 noise bias')
    pylab.errorbar(bins_centered,binstats_bfit_noisy['m2'],yerr=binstats_bfit_noisy['m2_std'], fmt='-bx', label = 'm2 noise bias')
    
    pylab.errorbar(bins_centered,binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'], fmt='-m+', label = 'm1 noise+model+interact')
    pylab.errorbar(bins_centered,binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'], fmt='-cx', label = 'm2 noise+model+interact')

    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    # pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)

    pylab.xlim([0,1.5])
    # pylab.ylim([-0.1, 0.1])

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))

    # pylab.yscale('symlog',linthreshy=0.05)
    pylab.xlabel('redshift zphot')
    pylab.ylabel('m_i')

    yticks = [-0.1,-0.05,-req2_m,-0.01,-req1_m,0.,req1_m,0.01,req2_m,0.05,0.1]
    pylab.yticks(yticks,[str(x) for x in  yticks ])
    pylab.ylim([-0.01,0.07])

    pylab.grid()
    pylab.legend(loc='upper left',ncol=2)


    filename_fig = 'figures/fig.bins.redshift.png'
    pylab.savefig(filename_fig)
    pylab.close()
    logger.info('saved %s' % filename_fig)



# size no cut ------------------------------------------------------------------------------------------------------------------


    binstats_real_noisy = tabletools.loadTable(filepath='bins.size_nocut.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.size_nocut.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.size_nocut.real.clear.cat',dtype=dtype_table_binstats)

    weights = binstats_real_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_real_noisy['bin_ngals']))
    logger.info( 'm1 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m1'] , weights ))
    logger.info( 'm2 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m2'] , weights ))

    weights = binstats_bfit_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_bfit_noisy['bin_ngals']))
    logger.info( 'm1 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m1'] , weights ))
    logger.info( 'm2 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m2'] , weights ))

    pylab.figure()
    pylab.clf()

    bins_centered = _binCentersSameLen(binstats_real_clear['bin_value'])

    pylab.errorbar(bins_centered,binstats_real_clear['m1'],yerr=binstats_real_clear['m1_std'], fmt='-y+', label = 'm1 model bias')
    pylab.errorbar(bins_centered,binstats_real_clear['m2'],yerr=binstats_real_clear['m2_std'], fmt='-g', label = 'm2 model bias')

    pylab.errorbar(bins_centered,binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'], fmt='-r+', label = 'm1 noise+model+interact')
    pylab.errorbar(bins_centered,binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'], fmt='-bx', label = 'm2 noise+model+interact')

    pylab.errorbar(bins_centered,binstats_bfit_noisy['m1'],yerr=binstats_bfit_noisy['m1_std'], fmt='-m+', label = 'm1 noise bias')
    pylab.errorbar(bins_centered,binstats_bfit_noisy['m2'],yerr=binstats_bfit_noisy['m2_std'], fmt='-cx', label = 'm2 noise bias')

    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    # print xadd
    # pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)
    # print pylab.xlim()

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))

    pylab.yscale('symlog',linthreshy=0.05)
    pylab.xlabel('size Rgp/Rp')
    pylab.ylabel('m_i')

    yticks = [-0.1,-0.05,-req2_m,-0.01,-req1_m,0.,req1_m,0.01,req2_m,0.05,0.1]
    pylab.yticks(yticks,[str(x) for x in  yticks ])

    pylab.grid()
    # pylab.legend(loc='lower right')


    pylab.xlim([1.,2.])
    # pylab.ylim([-1., 1.])

    filename_fig = 'figures/fig.bins.size_nocut.png'
    pylab.savefig(filename_fig)
    pylab.close()
    logger.info('saved %s' % filename_fig)

# size ------------------------------------------------------------------------------------------------------------------



    binstats_real_noisy = tabletools.loadTable(filepath='bins.size.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.size.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.size.real.clear.cat',dtype=dtype_table_binstats)

    weights = binstats_real_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_real_noisy['bin_ngals']))
    logger.info( 'm1 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m1'] , weights ))
    logger.info( 'm2 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m2'] , weights ))

    weights = binstats_bfit_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_bfit_noisy['bin_ngals']))
    logger.info( 'm1 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m1'] , weights ))
    logger.info( 'm2 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m2'] , weights ))

    pylab.figure()
    pylab.clf()

    bins_centered = _binCentersSameLen(binstats_real_clear['bin_value'])

    pylab.errorbar(bins_centered,binstats_real_clear['m1'],yerr=binstats_real_clear['m1_std'], fmt='-y+', label = 'm1 model bias')
    pylab.errorbar(bins_centered,binstats_real_clear['m2'],yerr=binstats_real_clear['m2_std'], fmt='-g', label = 'm2 model bias')

    pylab.errorbar(bins_centered,binstats_bfit_noisy['m1'],yerr=binstats_bfit_noisy['m1_std'], fmt='-m+', label = 'm1 noise bias')
    pylab.errorbar(bins_centered,binstats_bfit_noisy['m2'],yerr=binstats_bfit_noisy['m2_std'], fmt='-cx', label = 'm2 noise bias')

    pylab.errorbar(bins_centered,binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'], fmt='-r+', label = 'm1 noise+model+interact')
    pylab.errorbar(bins_centered,binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'], fmt='-bx', label = 'm2 noise+model+interact')

    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    # print xadd
    # pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)
    # print pylab.xlim()


    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))

    # pylab.yscale('symlog',linthreshy=0.05)
    pylab.xlabel('size Rgp/Rp')
    pylab.ylabel('m_i')

    yticks = [-0.1,-0.05,-req2_m,-0.01,-req1_m,0.,req1_m,0.01,req2_m,0.05,0.1]
    pylab.yticks(yticks,[str(x) for x in  yticks ])

    pylab.xlim([1.3,2.])
    pylab.ylim([-0.03, 0.06])

    pylab.grid()
    pylab.legend(loc='lower right',ncol=2)

    filename_fig = 'figures/fig.bins.size.png'
    pylab.savefig(filename_fig)
    pylab.close()
    logger.info('saved %s' % filename_fig)

# morphology ------------------------------------------------------------------------------------------------------------


    binstats_real_noisy = tabletools.loadTable(filepath='bins.modd.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.modd.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.modd.real.clear.cat',dtype=dtype_table_binstats)

    weights = binstats_real_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_real_noisy['bin_ngals']))
    logger.info( 'm1 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m1'] , weights ))
    logger.info( 'm2 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m2'] , weights ))

    weights = binstats_bfit_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_bfit_noisy['bin_ngals']))
    logger.info( 'm1 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m1'] , weights ))
    logger.info( 'm2 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m2'] , weights ))


    pylab.figure()
    pylab.clf()

    bins_centered = _binCentersSameLen(binstats_real_clear['bin_value'])


    pylab.errorbar(bins_centered,binstats_real_clear['m1'],yerr=binstats_real_clear['m1_std'], fmt='-y+', label = 'm1 model bias')
    pylab.errorbar(bins_centered,binstats_real_clear['m2'],yerr=binstats_real_clear['m2_std'], fmt='-g', label = 'm2 model bias')

    pylab.errorbar(bins_centered,binstats_bfit_noisy['m1'],yerr=binstats_bfit_noisy['m1_std'], fmt='-m+', label = 'm1 noise bias')
    pylab.errorbar(bins_centered,binstats_bfit_noisy['m2'],yerr=binstats_bfit_noisy['m2_std'], fmt='-cx', label = 'm2 noise bias')

    pylab.errorbar(bins_centered,binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'], fmt='-r+', label = 'm1 noise+model+interact')
    pylab.errorbar(bins_centered,binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'], fmt='-bx', label = 'm2 noise+model+interact')

    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    # pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)

    pylab.xlim([0.,27.5])
    # pylab.ylim([-0.1, 0.1])

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))

    pylab.yscale('symlog',linthreshy=0.05)
    pylab.xlabel('morph Hubble Seq')
    pylab.ylabel('m_i')

    yticks = [-0.15,-0.1,-0.05,-req2_m,-0.01,-req1_m,0.,req1_m,0.01,req2_m,0.05,0.1]
    pylab.yticks(yticks,[str(x) for x in  yticks ])
    pylab.ylim([-0.03,0.04])

    pylab.grid()
    pylab.legend(loc='lower right',ncol=2)


    filename_fig = 'figures/fig.bins.modd.png'
    pylab.savefig(filename_fig)
    pylab.close()
    logger.info('saved %s' % filename_fig)

# model bias ------------------------------------------------------------------------------------------------------------

    req1_m = 0.02
    req2_m = 0.004


    binstats_real_noisy = tabletools.loadTable(filepath='bins.model_bias.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.model_bias.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.model_bias.real.clear.cat',dtype=dtype_table_binstats)

    weights = binstats_real_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_real_noisy['bin_ngals']))
    logger.info( 'm1 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m1'] , weights ))
    logger.info( 'm2 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m2'] , weights ))

    weights = binstats_bfit_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_bfit_noisy['bin_ngals']))
    logger.info( 'm1 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m1'] , weights ))
    logger.info( 'm2 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m2'] , weights ))

    pylab.figure()
    pylab.clf()

    bins_centered = _binCentersSameLen(binstats_real_clear['bin_value'])

    pylab.errorbar(bins_centered,abs(binstats_real_clear['m1']),yerr=binstats_real_clear['m1_std'], fmt='-y+', label = 'm1 model bias')
    pylab.errorbar(bins_centered,abs(binstats_real_clear['m2']),yerr=binstats_real_clear['m2_std'], fmt='-g', label = 'm2 model bias')

    pylab.errorbar(bins_centered,binstats_bfit_noisy['m1'],yerr=binstats_bfit_noisy['m1_std'], fmt='-m+', label = 'm1 noise bias')
    pylab.errorbar(bins_centered,binstats_bfit_noisy['m2'],yerr=binstats_bfit_noisy['m2_std'], fmt='-cx', label = 'm2 noise bias')

    pylab.errorbar(bins_centered,binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'], fmt='-r+', label = 'm1 noise+model+interact')
    pylab.errorbar(bins_centered,binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'], fmt='-bx', label = 'm2 noise+model+interact')

    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    # pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)

    pylab.xlim([0.0,0.025])
    # pylab.ylim([-0.1, 0.1])

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))

    pylab.yscale('symlog',linthreshy=0.05)
    pylab.xlabel('m1')
    pylab.ylabel('m_i')

    yticks = [-0.15,-0.1,-0.05,-req2_m,-0.01,-req1_m,0.,req1_m,0.01,req2_m,0.05,0.1]
    pylab.yticks(yticks,[str(x) for x in  yticks ])
    pylab.ylim([-0.03,0.04])

    pylab.grid()
    pylab.legend(loc='lower right',ncol=2)


    filename_fig = 'figures/fig.bins.model_bias.png'
    pylab.savefig(filename_fig)
    pylab.close()
    logger.info('saved %s' % filename_fig)


def getCuts():

    cut_rgprp_min = 1.3
    cut_rgprp_max = 1000
    cut_snr_min = 40
    cut_snr_max = 100000

    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)

    psf_size = 0.765
    cut_stats = stats_array.copy()
    cut_ajs = ajs_array.copy()
    select = numpy.logical_and(cut_stats['rgp']/psf_size > cut_rgprp_min , cut_stats['rgp']/psf_size < cut_rgprp_max)
    logger.info('cutting on size using %7d galaxies out of %7d' % (sum(select),len(cut_stats)))
    cut_stats = cut_stats[select]
    cut_ajs = cut_ajs[select]
    select = numpy.logical_and(cut_stats['snr'] > cut_snr_min , cut_stats['snr'] < cut_snr_max)
    logger.info('cutting on snr  using %7d galaxies out of %7d' % (sum(select),len(cut_stats)))
    cut_stats = cut_stats[select]
    cut_ajs = cut_ajs[select]

    logger.info('total  cutting %7d galaxies out of %7d' % (len(cut_ajs),len(stats_array)))

    return cut_stats,cut_ajs

def testSaveBiasForBins():

    # load the data
    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)   
    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)

### ------------------------ galaxy size --- full (no cuts)

    # logging.info('galaxy size -- nocut')

    # psf_size = 0.7695
    # bin_column = stats_array['rgp']/psf_size
    # bins_size_test = [1,1.4,2]

    # do the real noisy case

    # logging.info('--------------------------- getting results for real noisy galaxies ---------------------------')
    # mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=stats_array['cosmos_id'],bin_values=bins_size_test,bin_column_name='test_bins_size')
    # filename_results = 'bins.test_bins_size.real.noisy.cat'
    # tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    # logging.info('--------------------------- getting results for bfit noisy galaxies ---------------------------')
    # mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_size_test,bin_column_name='test_bins_size')
    # filename_results = 'bins.test_bins_size.bfit.noisy.cat'
    # tabletools.saveTable(filename_results,mc_results)

    # # do the real noiseless case
    # logging.info('--------------------------- getting results for real galaxies ---------------------------')
    # mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_size_test,bin_column_name='test_bins_size')
    # filename_results = 'bins.test_bins_size.real.clear.cat'
    # tabletools.saveTable(filename_results,mc_results)



### ------------------------ galaxy model bias

    logging.info('galaxy model bias')

    cut_stats,cut_ajs = getCuts()
    stats_array = cut_stats;
    ajs_array   = cut_ajs;

    bin_column = abs(stats_array['m1']+1j*stats_array['m2'])
    bins_model_bias = [0.0, 0.005, 0.01, 0.015, 0.02 , 0.025] 

    
    # do the real noisy case
    logging.info('--------------------------- getting results for real noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=stats_array['cosmos_id'],bin_values=bins_model_bias,bin_column_name='bins_model_bias')
    filename_results = 'bins.model_bias.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    logging.info('--------------------------- getting results for bfit noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_model_bias,bin_column_name='bins_model_bias')
    filename_results = 'bins.model_bias.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the real noiseless case
    logging.info('--------------------------- getting results for real galaxies ---------------------------')
    mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_model_bias,bin_column_name='bins_model_bias')
    filename_results = 'bins.model_bias.real.clear.cat'
    tabletools.saveTable(filename_results,mc_results)

# model bias ------------------------------------------------------------------------------------------------------------

    req1_m = 0.02
    req2_m = 0.004


    binstats_real_noisy = tabletools.loadTable(filepath='bins.model_bias.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.model_bias.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.model_bias.real.clear.cat',dtype=dtype_table_binstats)

    weights = binstats_real_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_real_noisy['bin_ngals']))
    logger.info( 'm1 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m1'] , weights ))
    logger.info( 'm2 real noisy sum check %2.5f' % numpy.inner(binstats_real_noisy['m2'] , weights ))

    weights = binstats_bfit_noisy['bin_ngals'].astype(numpy.float64) / float(sum(binstats_bfit_noisy['bin_ngals']))
    logger.info( 'm1 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m1'] , weights ))
    logger.info( 'm2 bfit noisy sum check %2.5f' % numpy.inner(binstats_bfit_noisy['m2'] , weights ))

    pylab.figure()
    pylab.clf()

    bins_centered = _binCentersSameLen(binstats_real_clear['bin_value'])

    pylab.errorbar(bins_centered,abs(binstats_real_clear['m1']),yerr=binstats_real_clear['m1_std'], fmt='-y+', label = 'm1 model bias')
    pylab.errorbar(bins_centered,abs(binstats_real_clear['m2']),yerr=binstats_real_clear['m2_std'], fmt='-g', label = 'm2 model bias')

    pylab.errorbar(bins_centered,binstats_bfit_noisy['m1'],yerr=binstats_bfit_noisy['m1_std'], fmt='-m+', label = 'm1 noise bias')
    pylab.errorbar(bins_centered,binstats_bfit_noisy['m2'],yerr=binstats_bfit_noisy['m2_std'], fmt='-cx', label = 'm2 noise bias')

    pylab.errorbar(bins_centered,binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'], fmt='-r+', label = 'm1 noise+model+interact')
    pylab.errorbar(bins_centered,binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'], fmt='-bx', label = 'm2 noise+model+interact')

    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    # pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)

    pylab.xlim([-0.005,0.025])
    # pylab.ylim([-0.1, 0.1])

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))

    pylab.yscale('symlog',linthreshy=0.05)
    pylab.xlabel('m1')
    pylab.ylabel('m_i')

    yticks = [-0.15,-0.1,-0.05,-req2_m,-0.01,-req1_m,0.,req1_m,0.01,req2_m,0.05,0.1]
    pylab.yticks(yticks,[str(x) for x in  yticks ])
    pylab.ylim([-0.03,0.04])

    pylab.grid()
    pylab.legend(loc='lower right',ncol=2)


    filename_fig = 'figures/fig.bins.model_bias.png'
    pylab.savefig(filename_fig)
    pylab.close()




def saveBiasForBins():

    # load the data
    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)   
    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)



### ------------------------ galaxy size --- full (no cuts)

    logging.info('galaxy size -- nocut')


    psf_size = 0.7695
    bin_column = stats_array['rgp']/psf_size

    # do the real noisy case

    logging.info('--------------------------- getting results for real noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=stats_array['cosmos_id'],bin_values=bins_size_nocut,bin_column_name='size_nocut_real_noisy')
    filename_results = 'bins.size_nocut.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    logging.info('--------------------------- getting results for bfit noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_size_nocut,bin_column_name='size_nocut_bfit_noisy')
    filename_results = 'bins.size_nocut.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # # do the real noiseless case
    logging.info('--------------------------- getting results for real galaxies ---------------------------')
    mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_size_nocut,bin_column_name='size_nocut_real_clear')
    filename_results = 'bins.size_nocut.real.clear.cat'
    tabletools.saveTable(filename_results,mc_results)


### ------------------------ size cut

    cut_stats,cut_ajs = getCuts()
    stats_array = cut_stats;
    ajs_array   = cut_ajs;

### ------------------------ redshift

    logging.info('galaxy redshift')

    # bins_redshift = [0 , 0.35, 0.6, 0.8, 1.1 , 1.5]    

    # do the real noisy case
    logging.info('--------------------------- getting results for real noisy galaxies ---------------------------')
    mc_results_real_noisy = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=ajs_array['ZPHOT'], bin_ids=ajs_array['IDENT'],bin_values=bins_redshift,bin_column_name='zphot_real_noisy')
    filename_results = 'bins.redshift.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results_real_noisy)

    # do the bfit noisy case
    logging.info('--------------------------- getting results for bfit noisy galaxies ---------------------------')
    mc_results_bfit_noisy = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=ajs_array['ZPHOT'], bin_ids=ajs_array['IDENT'],bin_values=bins_redshift,bin_column_name='zphot_bfit_noisy')
    filename_results = 'bins.redshift.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results_bfit_noisy)

    # do the real noiseless case
    logging.info('--------------------------- getting results for real galaxies ---------------------------')
    mc_results_real_clear = getBiasForBins(results_real,truth_array_26000,bin_column=ajs_array['ZPHOT'], bin_ids=ajs_array['IDENT'],bin_values=bins_redshift,bin_column_name='zphot_bfit_clear')
    filename_results = 'bins.redshift.real.clear.cat'
    tabletools.saveTable(filename_results,mc_results_real_clear)


### ------------------------ galaxy size -- cut

    logging.info('galaxy size -- cut')


    psf_size = 0.7695
    bin_column = stats_array['rgp']/psf_size

    # do the real noisy case

    logging.info('--------------------------- getting results for real noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=stats_array['cosmos_id'],bin_values=bins_size,bin_column_name='size_real_noisy')
    filename_results = 'bins.size.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    logging.info('--------------------------- getting results for bfit noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_size,bin_column_name='size_bfit_noisy')
    filename_results = 'bins.size.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # # do the real noiseless case
    logging.info('--------------------------- getting results for real galaxies ---------------------------')
    mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_size,bin_column_name='size_real_clear')
    filename_results = 'bins.size.real.clear.cat'
    tabletools.saveTable(filename_results,mc_results)




### ------------------------ galaxy morphology


    logging.info('galaxy morphology')


    # bins_modd = [0 , 5 ,  9 , 20, 32]
    bin_column = ajs_array['MODD']
    bin_ids = stats_array['cosmos_id']

        # do the real noisy case

    logging.info('--------------------------- getting results for %d real noisy galaxies ---------------------------' % len(bin_ids))
    mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=bin_ids,bin_values=bins_modd,bin_column_name='modd_real_noisy')
    filename_results = 'bins.modd.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    logging.info('--------------------------- getting results for %d bfit noisy galaxies ---------------------------'  % len(bin_ids))
    mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=bin_ids,bin_values=bins_modd,bin_column_name='modd_bfit_noisy')
    filename_results = 'bins.modd.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # # do the real noiseless case
    logging.info('--------------------------- getting results for %d real galaxies ---------------------------'  % len(bin_ids))
    mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=bin_ids,bin_values=bins_modd,bin_column_name='modd_real_clear')
    filename_results = 'bins.modd.real.clear.cat'
    tabletools.saveTable(filename_results,mc_results)



### ------------------------ galaxy model bias

    logging.info('galaxy model bias')

    cut_stats,cut_ajs = getCuts()
    stats_array = cut_stats;
    ajs_array   = cut_ajs;

    bin_column = abs(stats_array['m1']+1j*stats_array['m2'])
    bins_model_bias = [0.0, 0.005, 0.01, 0.015, 0.02 , 0.025] 

    
    # do the real noisy case
    logging.info('--------------------------- getting results for real noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=stats_array['cosmos_id'],bin_values=bins_model_bias,bin_column_name='bins_model_bias')
    filename_results = 'bins.model_bias.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    logging.info('--------------------------- getting results for bfit noisy galaxies ---------------------------')
    mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_model_bias,bin_column_name='bins_model_bias')
    filename_results = 'bins.model_bias.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the real noiseless case
    logging.info('--------------------------- getting results for real galaxies ---------------------------')
    mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=stats_array['cosmos_id'],bin_values=bins_model_bias,bin_column_name='bins_model_bias')
    filename_results = 'bins.model_bias.real.clear.cat'
    tabletools.saveTable(filename_results,mc_results)





def _binCenters(b):
      c1 = [ b[i-1] + (b[i] - b[i-1])/2. for i in range(1,len(b)) ] 
      return c1

def _binCentersSameLen(b):

    c1 = [ b[i] + (b[i+1] - b[i])/2. for i in range(0,len(b)-1) ] 
    c1.append(b[-1] + (b[-1] - b[-2])/2.)
    return c1


def plotModelBias():

    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    results_stats      = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)
    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)
    
    results_modd = ajs_array['MODD']

    pixel_scale = 0.27
    n_bins = 20
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
    results_modd = results_modd[select]

    # pylab.ion()
    pylab.close('all')

    # m1 m2 hist    

    pylab.figure(figsize=(5,4))
    pylab.clf()
    n_bins = 50
    h1,b1,_ = pylab.hist(results_m1,n_bins,color='r',histtype='step',label='m1')
    h2,b2,_ = pylab.hist(results_m2,n_bins,color='b',histtype='step',label='m2')
    # pylab.plot(_binCenters(b1),h1,'r+-')
    # pylab.plot(_binCenters(b2),h2,'bx-')
    pylab.axvline(x=results_m1.mean(),linewidth=1, color='r')
    pylab.axvline(x=results_m2.mean(),linewidth=1, color='b')

    pylab.xlabel('multiplicative bias')
    pylab.grid()
    pylab.legend()
    print "mean" , numpy.mean(results_m1)
    print "mean" , numpy.mean(results_m2)
    print "stdv" , numpy.std(results_m1,ddof=1)
    print "stdv" , numpy.std(results_m2,ddof=1)
    print "stdm" , numpy.std(results_m1,ddof=1)/numpy.sqrt(len(results_m1))
    print "stdm" , numpy.std(results_m2,ddof=1)/numpy.sqrt(len(results_m2))
    print 'nans?' , any(numpy.isnan(results_m1))
    filename_fig = 'figures/figure.hist.m1m2.real.png'
    pylab.savefig(filename_fig)
    filename_fig = 'figures/figure.hist.m1m2.real.eps'
    pylab.savefig(filename_fig,dpi=1000)
    logger.info('saved %s' % filename_fig)

    # snr distribution

    pylab.figure()
    h,b,p = pylab.hist(results_snr,bins=bins_snr)
    pylab.xlabel('snr')
    pylab.xscale('log')    
    filename_fig = 'figures/figure.hist.snr.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # Rgp/Rp distribution
    
    pylab.figure()
    psf_size =  0.7695
    results_size = results_rgp/psf_size
    h,b,p = pylab.hist(results_size,bins=bins_size)
    pylab.xlabel('rgp/rp')    
    filename_fig = 'figures/figure.hist.size.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # hlr distribution

    pylab.figure()
    results_hlr = results_stats['hlr']/pixel_scale
    h,b,p = pylab.hist(results_hlr,bins=bins_hlr)
    pylab.xlabel('hlr [pixels]')    
    filename_fig = 'figures/figure.hist.hlr.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # redshift distribution

    pylab.figure()
    h,b,p = pylab.hist(results_zphot,bins=n_bins)
    # h,b,p = pylab.hist(results_zphot,bins=bins_redshift)
    pylab.xlabel('photo-z')    
    filename_fig = 'figures/figure.hist.redshift.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # modd distribution

    pylab.figure()
    # bins_modd = [s for s in set(results_modd)]
    # bins_modd = [0,9,32]
    h,b,p = pylab.hist(results_modd,bins=bins_modd)
    # n_bins_modd = 31
    # h,b,p = pylab.hist(results_modd,bins=n_bins)
    pylab.xlabel('modd')    
    filename_fig = 'figures/figure.hist.modd.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)


    # m1 vs snr

    digitized = numpy.digitize(results_snr, bins_snr)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(0, len(bins_snr))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(0, len(bins_snr))]
    bin_stdm1 = [results_m1[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m1[digitized == i])) for i in range(0, len(bins_snr))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(0, len(bins_snr))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(0, len(bins_snr))]
    bin_stdm2 = [results_m1[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m2[digitized == i])) for i in range(0, len(bins_snr))]

    
    pylab.figure()
    # pylab.errorbar(bins_snr,bin_mean1,yerr=bin_stdv1,fmt='r')
    # pylab.errorbar(bins_snr,bin_mean2,yerr=bin_stdv2,fmt='b')
    pylab.errorbar(bins_snr,bin_mean1,yerr=bin_stdm1,fmt='m')
    pylab.errorbar(bins_snr,bin_mean2,yerr=bin_stdm2,fmt='c')
    pylab.xscale('log')
    pylab.xlabel('snr')
    pylab.ylabel('m1')
    filename_fig = 'figures/figure.hist.M1M2vsSNR.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # m1 vs size


    digitized = numpy.digitize(results_size, bins_size)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(0, len(bins_size))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(0, len(bins_size))]
    bin_stdm1 = [results_m1[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m1[digitized == i])) for i in range(0, len(bins_size))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(0, len(bins_size))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(0, len(bins_size))]
    bin_stdm2 = [results_m2[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m2[digitized == i])) for i in range(0, len(bins_size))]
    
    pylab.figure()
    # pylab.errorbar(bins_size,bin_mean1,yerr=bin_stdv1,fmt='r')
    # pylab.errorbar(bins_size,bin_mean2,yerr=bin_stdv2,fmt='b')
    pylab.errorbar(bins_size,bin_mean1,yerr=bin_stdm1,fmt='m')
    pylab.errorbar(bins_size,bin_mean2,yerr=bin_stdm2,fmt='c')
    pylab.xlabel('size')
    pylab.ylabel('m1')

    filename_fig = 'figures/figure.hist.M1M2vsSIZE.real.png'
    pylab.savefig(filename_fig)

    # snr vs redshifts

    pylab.figure() 
    digitized = numpy.digitize(results_zphot, bins_redshift)
    for i in range(1, len(bins_redshift)):
        h,b = pylab.histogram(results_snr[digitized == i],bins=bins_snr,normed=True)
        pylab.plot(_binCenters(b),h,'x-',label='redshift bin %2.2f' % bins_redshift[i])

    pylab.legend()
    pylab.xlabel('snr')
    pylab.ylabel('normalised hist')


    filename_fig = 'figures/figure.hist.SNRvsREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # size vs redshifts

    pylab.figure() 
    bins_size_extra = numpy.linspace(1,4,20)
    digitized = numpy.digitize(results_zphot, bins_redshift)
    for i in range(1, len(bins_redshift)):
        h,b = pylab.histogram(results_size[digitized == i],bins=bins_size_extra,normed=True)
        # pylab.hist(results_size[digitized == i],bins=bins_size_extra,normed=True,histtype='step')
        pylab.plot(_binCenters(b),h,label='redshift bin %2.2f' % bins_redshift[i])
    pylab.legend()
    pylab.xlabel('Rgp/Rp')
    pylab.ylabel('normalised hist')

    filename_fig = 'figures/figure.hist.SIZEvsREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # redshifts vs redshifts

    pylab.figure() 
    digitized = numpy.digitize(results_zphot, bins_redshift)
    for i in range(1, len(bins_redshift)):
        h,b = pylab.histogram(results_zphot[digitized == i],bins=bins_redshift,normed=True)
        pylab.plot(_binCenters(b),h,'x',label='redshift bin %2.2f' % bins_redshift[i])
    pylab.legend()
    pylab.xlabel('zphot - control')
    pylab.ylabel('normalised hist')

    # m1 vs modd

    digitized = numpy.digitize(results_modd , bins_modd)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(1, len(bins_modd))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(1, len(bins_modd))]
    bin_stdm1 = [results_m1[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m1[digitized == i])) for i in range(1, len(bins_modd))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(1, len(bins_modd))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(1, len(bins_modd))]
    bin_stdm2 = [results_m2[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m2[digitized == i])) for i in range(1, len(bins_modd))]
    
    pylab.figure()
    pylab.errorbar(bins_modd[1:],bin_mean1,yerr=bin_stdv1)
    pylab.errorbar(bins_modd[1:],bin_mean2,yerr=bin_stdv2)
    pylab.errorbar(bins_modd[1:],bin_mean1,yerr=bin_stdm1)
    pylab.errorbar(bins_modd[1:],bin_mean2,yerr=bin_stdm2)
    pylab.xlabel('modd')
    pylab.ylabel('m1')
    filename_fig = 'figures/figure.hist.M1M2vsMODD.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    pylab.figure() 
    digitized = numpy.digitize(results_modd, bins_modd)
    n_bins_m = 100


    colors_modd_m1 = getColorMap(len(bins_modd))
    colors_modd_m2 = getColorMap(len(bins_modd))

    # colors_modd_m1 = ['r','b']
    # colors_modd_m2 = ['m','c']
    for i in range(1, len(bins_modd)):
        m1= results_m1[digitized == i]
        m2= results_m2[digitized == i]

        h,b,p = pylab.hist(m1,bins=n_bins_m,normed=True,histtype='step',label='m1 Hubble Seq %d-%d' % (bins_modd[i-1],bins_modd[i]),color=colors_modd_m1[i-1])
        h,b,p = pylab.hist(m2,bins=n_bins_m,normed=True,histtype='step',label='m2 Hubble Seq %d-%d' % (bins_modd[i-1],bins_modd[i]),color=colors_modd_m2[i-1])
        pylab.axvline(linewidth=2, color='k')
        pylab.axvline(x=m1.mean(),linewidth=1, color=colors_modd_m1[i-1])
        pylab.axvline(x=m2.mean(),linewidth=1, color=colors_modd_m2[i-1])

    pylab.grid()
    pylab.legend()
    pylab.xlabel('m1m2')
    pylab.xticks([-0.1 , -0.05 , -0.01 ,   0.01 , 0.05 , 0.1])
    pylab.ylabel('normalised hist')

    filename_fig = 'figures/figure.hist.M1M2forMODD.real.eps'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # m1 vs redshift

    digitized = numpy.digitize(results_zphot, bins_redshift)
    bin_mean1 = [results_m1[digitized == i].mean()      for i in range(1, len(bins_redshift))]
    bin_stdv1 = [results_m1[digitized == i].std(ddof=1) for i in range(1, len(bins_redshift))]
    bin_stdm1 = [results_m1[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m1[digitized == i])) for i in range(1, len(bins_redshift))]
    bin_mean2 = [results_m2[digitized == i].mean()      for i in range(1, len(bins_redshift))]
    bin_stdv2 = [results_m2[digitized == i].std(ddof=1) for i in range(1, len(bins_redshift))]
    bin_stdm2 = [results_m2[digitized == i].std(ddof=1)/numpy.sqrt(len(results_m2[digitized == i])) for i in range(1, len(bins_redshift))]
    
    pylab.figure()
    # pylab.errorbar(bins_redshift[1:],bin_mean1,yerr=bin_stdv1)
    # pylab.errorbar(bins_redshift[1:],bin_mean2,yerr=bin_stdv2)
    pylab.errorbar(bins_redshift[1:],bin_mean1,yerr=bin_stdm1)
    pylab.errorbar(bins_redshift[1:],bin_mean2,yerr=bin_stdm2)
    pylab.xlabel('redshift')
    pylab.ylabel('m1')
    
    filename_fig = 'figures/figure.hist.M1M2vsREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    print len(results_stats)
    print len(results_size[results_size > 1.2])

    # redshift cross - sections

    pylab.figure()
    n_bins_m = 32

    redshift_color_tuples = getColorMap(len(bins_redshift))

    for i in range(1, len(bins_redshift)):
        m1= results_m1[digitized == i]
        m2= results_m2[digitized == i]
        h,b,p = pylab.hist(m1,bins=n_bins_m,normed=True,histtype='step',label='m1,2 z=%2.2f' % (bins_redshift[i]),color=redshift_color_tuples[i-1])
        h,b,p = pylab.hist(m2,bins=n_bins_m,normed=True,histtype='step', color=redshift_color_tuples[i-1])
        pylab.axvline(linewidth=2, color='k')
        pylab.axvline(x=m1.mean(),linewidth=1, color=redshift_color_tuples[i-1])
        pylab.axvline(x=m2.mean(),linewidth=1, color=redshift_color_tuples[i-1])        

    pylab.grid()
    pylab.legend()
    pylab.xlabel('m1m2')
    pylab.xticks([-0.1 , -0.05 , -0.01 ,   0.01 , 0.05 , 0.1])
    pylab.ylabel('normalised hist')

    filename_fig = 'figures/figure.hist.M1M2forREDSHIFT.real.png'
    pylab.savefig(filename_fig,dpi=200)
    logger.info('saved %s' % filename_fig)

    # modd vs redshift

    pylab.figure() 
    digitized = numpy.digitize(results_modd, bins_modd)
    for i in range(1, len(bins_modd)):
        # h,b,p = pylab.hist(results_zphot[digitized == i],bins=bins_redshift,normed=True,histtype='step',label='type=%d' % i )
        h,b = pylab.histogram(results_zphot[digitized == i],bins=bins_redshift,normed=True)
        pylab.plot(_binCenters(b),h,'-x',label='Hubble Seq=%d-%d' % (bins_modd[i-1],bins_modd[i]))

        
    pylab.grid()
    pylab.legend()
    pylab.xlabel('redshift')
    pylab.ylabel('normalised hist')

    filename_fig = 'figures/figure.hist.MMODforREDSHIFT.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # modd vs size

    pylab.figure() 
    digitized = numpy.digitize(results_modd, bins_modd)
    bins_size_extra = numpy.arange(1 , 1.8 , 0.05)
    for i in range( 1, len(bins_modd) ):
        # print set(digitized)
        # h,b,p = pylab.hist(results_size[digitized == i],bins=bins_size_extra,normed=True,histtype='step',label='Hubble Seq=%d-%d' % (bins_modd[i-1],bins_modd[i]))
        # h,b,p = pylab.hist(results_size[digitized == i],bins=bins_size_extra,normed=True,histtype='bar',label='type=%d' % i , alpha = 0.1)
        h,b = pylab.histogram(results_size[digitized == i],bins=bins_size_extra,normed=True)
        pylab.plot(_binCenters(b),h,'-x',label='Hubble Seq=%d-%d' % (bins_modd[i-1],bins_modd[i]))
        
    pylab.grid()
    pylab.legend()
    pylab.xlabel('size Rgp/Rp')
    pylab.ylabel('normalised hist')
    pylab.xlim([1,1.8])

    filename_fig = 'figures/figure.hist.MMODforSIZE.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

def plotTime():

    # load the data
    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)   
    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)


    results_bfit_noisy_valid = results_bfit_noisy[results_bfit_noisy['time_taken']<666]
    results_real_noisy_valid = results_real_noisy[results_real_noisy['time_taken']<666]

    # split = len(results_bfit_noisy_valid)/2

    # pylab.figure()
    # # pylab.plot(results_bfit_noisy_valid['time_taken'],'x')
    # # pylab.show()
    # pylab.hist(results_bfit_noisy_valid['time_taken'][0:split],1000,histtype='step',color='r')
    # pylab.hist(results_bfit_noisy_valid['time_taken'][split:] ,1000,histtype='step',color='b')
    # pylab.show()


    # pylab.figure()
    # # pylab.plot(results_real_noisy_valid['time_taken'],'x')
    # # pylab.show()
    # pylab.hist(results_real_noisy_valid['time_taken'][0:split],1000,histtype='step',color='r')
    # pylab.hist(results_real_noisy_valid['time_taken'][split:] ,1000,histtype='step',color='b')
    # pylab.show()

    n_per_task = 640
    n_tasks = len(results_bfit_noisy)/n_per_task
    time_mean = []
    time_total = []
    for i in range(0,len(results_bfit_noisy),n_per_task):
        i_start = i
        i_end   = i + n_per_task
        time_mean.append(numpy.mean(results_bfit_noisy['time_taken'][i_start:i_end]))
        time_total.append(numpy.sum(results_bfit_noisy['time_taken'][i_start:i_end]))

    time_mean  = numpy.array(time_mean)
    time_total = numpy.array(time_total)
    time_mean_valid  = time_mean[time_mean < NO_RESULT_FLAG]  
    time_total_valid = time_total[time_mean < NO_RESULT_FLAG]

    pylab.hist(time_mean_valid,100)
    pylab.xlabel('mean time per galaxy [sec]')
    pylab.show()
    pylab.hist(time_total_valid/60.,100)
    pylab.xlabel('time total per file [min]')
    pylab.show()








def getBiasForBins(results_array,truth_array,bin_column, bin_ids,bin_values , bin_column_name='test'):

    """
    @brief plot bias for different bins of different parameters
    @results_array use this array to get the results -- it can be noiseless, noisy real or noisy bfit
    @truth_array truth array used to create results_array
    @info_string some info for plot titles
    """

    # initialise bins
    bins_ids_split = []
    digitized = numpy.digitize(bin_column, bin_values)
    bins_ids_split = [bin_ids[digitized == i] for i in range(1, len(bin_values))] 
    
    mc_results = []

    n_gals_in_all_bins = 0


    # loop over bins
    for i,ids in enumerate(bins_ids_split):
        results_bin,truth_bin,_ = analyse.selectByIDs(ids,results_array,truth_array,logger)
        n_gals_in_all_bins += len(results_bin)
        logger.info('bin %d, number of galaxies in sample %5d, n_gals_in_all_bins=%d' % (i,len(results_bin),n_gals_in_all_bins))
        mcr = analyse.getBiasForResults(results_bin,truth_bin,logger=logger,bin_param=bin_column_name,bin_id=i, n_gals_per_mean=10000)
        mcr['bin_id'] = i
        mcr['bin_value'] = bin_values[i]
        mc_results.append(mcr)

        if logger.level <= logging.DEBUG:
            import pylab

            results_bin_use = results_bin[results_bin['e1'] < 1]

            pylab.figure()
            pylab.hist(results_bin_use['e1'],color='r',bins=200,histtype='step')
            pylab.hist(results_bin_use['e2'],color='b',bins=200,histtype='step')
            pylab.xlabel('e1,e2')
            pylab.ylabel('count')
            pylab.title('bins pe %s %02d %f ' % (bin_column_name,i,bin_values[i]))
            filename_fig = './debug/fig.bins.pe.%s.%02d.png' % (bin_column_name,i)
            pylab.savefig(filename_fig)
            pylab.close()


    mc_results = numpy.concatenate(mc_results)

    return mc_results

   
def main():

    global logger , config , args

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('command', type=str, help='what to do?')
    # parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    # parser.add_argument('--filepath_truth', type=str, default='truth.26000.sorted.pp', help='truth file for the run, overrides the config file (by default is taken from yaml file)')
    # parser.add_argument('--filepath_stats', type=str, default='stats.nmb_main.real.pp', help='stats file')
    # parser.add_argument('--filepath_results', type=str, default='results.nmb_main.real.pp', help='results file')
    # parser.add_argument('--filepath_acs', type=str, default='cosmos_acs_shera_may2011.fits.gz', help='cosmos_acs_shera_may2011.fits.gz file')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')

    args = parser.parse_args()
    # args.name_config = os.path.basename(args.filepath_config).replace('.yaml','')

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    logger = logging.getLogger("nmb_main_plots.py") 
    logger.setLevel(logging_level)

    if not os.path.isdir('figures'): os.mkdir('figures')
    if (args.verbosity > 2) and (not os.path.isdir('figures')): os.mkdir('figures')
    
    # load the configuration file
    # config = yaml.load(open(args.filepath_config,'r')) 
    # store the args in config so it's easier to use them
    # config['args'] = args    

    # plotModelBias()
    # getTableACSjoinStats()

    # get total bias

    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)
    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)
   
    # global results_bfit_noisy
    if results_bfit_noisy['id_cosmos'][0] == 1000270:
        logger.warning('fixing bugs with ID')
        results_bfit_noisy_bug = tabletools.loadTable( filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger); 
        results_bfit_noisy = results_bfit_noisy_bug.copy()
        results_bfit_noisy['id_cosmos'] = results_bfit_noisy['id_cosmos'] / 10
        tabletools.addTable(table_name='results_bfit_noisy',table=results_bfit_noisy)

    # cut_stats,cut_ajs = getCuts()
    # stats_array = cut_stats;
    # ajs_array   = cut_ajs;
    # ids = ajs_array['IDENT']

    # results_bfit_noisy_cut,truth_bfit_cut,_ = analyse.selectByIDs(ids,results_bfit_noisy,truth_array_25880,logger)
    # results_real_noisy_cut,truth_real_cut,_ = analyse.selectByIDs(ids,results_real_noisy,truth_array_25880,logger)

    # mc_bfit_noisy = analyse.getBiasForResults(results_bfit_noisy_cut,truth_bfit_cut,n_gals_per_mean=2e6,bin_param='bfit_all_cut',bin_id=0,logger=logger)
    # mc_real_noisy = analyse.getBiasForResults(results_real_noisy_cut,truth_real_cut,n_gals_per_mean=2e6,bin_param='real_all_cut',bin_id=0,logger=logger)


    # print 'bfit' , mc_bfit_noisy
    # print 'real' , mc_real_noisy

    # mc_bfit_noisy = analyse.getBiasForResults(results_bfit_noisy,truth_array_25880,n_gals_per_mean=2e6,bin_param='bfit_all',bin_id=0,logger=logger)
    # mc_real_noisy = analyse.getBiasForResults(results_real_noisy,truth_array_25880,n_gals_per_mean=2e6,bin_param='real_all',bin_id=0,logger=logger)

    # print 'bfit' , mc_bfit_noisy
    # print 'real' , mc_real_noisy

    # get the bins
    
    # testSaveBiasForBins()
    # saveBiasForBins()
    # plotBiasForBins()
    # plotTime()

    eval(args.command + '()')

    

if __name__ == "__main__":

    main()