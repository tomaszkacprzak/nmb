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
filepath_truth_26000        = 'truth.26000.cat'
filepath_acs                = 'cosmos_acs_shera_may2011.fits.gz'
filepath_results_real       = 'results.nmb_main.real.pp' 
filepath_results_real_noisy = 'results.nmb_main.real.noisy.fits'
filepath_results_bfit_noisy = 'results.nmb_main.bfit.noisy.fits'

NO_RESULT_FLAG = 666

def getTableACSjoinStats():

    global results_stats

    file_pickle = open(args.filepath_stats)
    results_stats = pickle.load(file_pickle)
    logger.info('opened stats file %s with %d rows' % (args.filepath_stats,results_stats.shape[0]))

    logger.info('getting modd')
    acs_array = pyfits.getdata(args.filepath_acs,1)
    results_modd = numpy.zeros(len(results_stats))
    for i in range(len(results_modd)):
        select = numpy.nonzero(acs_array['IDENT']==results_stats['cosmos_id'][i])
        results_modd[i] = acs_array['modd'][select]
    logger.info('got modd')
    tabletools.saveTable(filename_acs_join_stats,results_modd)
    return results_modd

def plotBiasForBins():

    req1_m = 0.02
    req2_m = 0.003

    # redshift

    binstats_real_noisy = tabletools.loadTable(filepath='bins.redshift.real.noisy.cat',dtype=dtype_table_binstats)
    binstats_bfit_noisy = tabletools.loadTable(filepath='bins.redshift.bfit.noisy.cat',dtype=dtype_table_binstats)
    binstats_real_clear = tabletools.loadTable(filepath='bins.redshift.real.clear.cat',dtype=dtype_table_binstats)

    pylab.figure()
    pylab.clf()

    pylab.errorbar(binstats_real_clear['bin_value'],binstats_real_clear['m1'],yerr=binstats_real_clear['m1_std'])
    pylab.errorbar(binstats_real_clear['bin_value'],binstats_real_clear['m2'],yerr=binstats_real_clear['m2_std'])

    pylab.errorbar(binstats_real_noisy['bin_value'],binstats_real_noisy['m1'],yerr=binstats_real_noisy['m1_std'])
    pylab.errorbar(binstats_real_noisy['bin_value'],binstats_real_noisy['m2'],yerr=binstats_real_noisy['m2_std'])
    xadd = max( [ abs(binstats_real_noisy['bin_value'].min()) , abs(binstats_real_noisy['bin_value'].max()) ] ) *0.1
    print xadd
    pylab.xlim(binstats_real_noisy['bin_value'].min()-xadd,binstats_real_noisy['bin_value'].max()+xadd)
    print pylab.xlim()

    corner = pylab.xlim()[0]
    length = abs(pylab.xlim()[1]) + abs(pylab.xlim()[0])
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req1_m), length , 2*req1_m , alpha=0.1))
    pylab.gca().add_patch(pylab.Rectangle(  (corner,-req2_m), length , 2*req2_m , alpha=0.2))


    pylab.xlabel('redshift zphot')
    pylab.ylabel('m_i')


    filename_fig = 'figures/fig.bins.redshift.png'
    pylab.savefig(filename_fig)
    pylab.close()








def saveBiasForBins():

    # load the data
    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)
    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)

    # m vs redshift


    ### ------------------------ redshift

    # bins_redshift = [0 , 0.35, 0.6, 0.8, 1.1 , 1.5]    

    # # do the real noisy case
    # mc_results_real_noisy = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=acs_array['zphot'], bin_ids=acs_array['IDENT'],bin_values=bins_redshift,bin_column_name='zphot_real_noisy')
    # filename_results = 'bins.redshift.real.noisy.cat'
    # tabletools.saveTable(filename_results,mc_results_real_noisy)

    # # do the bfit noisy case
    # mc_results_bfit_noisy = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=acs_array['zphot'], bin_ids=acs_array['IDENT'],bin_values=bins_redshift,bin_column_name='zphot_bfit_noisy')
    # filename_results = 'bins.redshift.bfit.noisy.cat'
    # tabletools.saveTable(filename_results,mc_results_bfit_noisy)

    # # do the real noiseless case
    # mc_results_real_clear = getBiasForBins(results_real,truth_array_26000,bin_column=acs_array['zphot'], bin_ids=acs_array['IDENT'],bin_values=bins_redshift,bin_column_name='zphot_bfit_clear')
    # filename_results = 'bins.redshift.real.clear.cat'
    # tabletools.saveTable(filename_results,mc_results_real_clear)

    ### ------------------------ galaxy size

    # n_bins_size = 12;
    # bins_size = numpy.linspace(0,3,n_bins_size)
    bins_size = [1 , 1.4 , 2.6]
    psf_size = 0.7695
    bin_column = stats_array['rgp']/psf_size
    import pdb;pdb.set_trace()
    # import pylab;pylab.hist(bin_column,bins_size);pylab.show()

    # do the real noisy case

    mc_results = getBiasForBins(results_real_noisy,truth_array_25880,bin_column=bin_column,bin_ids=truth_array_25880['id_cosmos'],bin_values=bins_size,bin_column_name='size_real_noisy')
    filename_results = 'bins.size.real.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the bfit noisy case
    mc_results = getBiasForBins(results_bfit_noisy,truth_array_25880,bin_column=bin_column, bin_ids=truth_array_25880['id_cosmos'],bin_values=bins_size,bin_column_name='size_bfit_noisy')
    filename_results = 'bins.size.bfit.noisy.cat'
    tabletools.saveTable(filename_results,mc_results)

    # do the real noiseless case
    mc_results = getBiasForBins(results_real,truth_array_26000,bin_column=bin_column, bin_ids=truth_array_25880['id_cosmos'],bin_values=bins_size,bin_column_name='size_real_clear')
    filename_results = 'bins.size.real.cat'
    tabletools.saveTable(filename_results,mc_results)


    # m vs size

    # n_bins_size = 5
    # bins_size = numpy.linspace(1,4,n_bins_size)
    # mc_results = getBiasForBins(results_real_noisy,truth_array,bin_column=stats_array['rgp'], bin_ids=stats_array['id_cosmos'],bin_values=bins_redshift)
    # filename_results = 'bins.redshift.cat'
    # writeMCresults(filename_results,mc_results)
    # # do the bfit noisy case
    # mc_results = getBiasForBins(results_bfit_noisy,truth_array,bin_column=stats_array['rgp'], bin_ids=stats_array['id_cosmos'],bin_values=bins_redshift)
    # filename_results = 'bins.redshift.bfit.noisy.cat'
    # writeMCresults(filename_results,mc_results)
    # # do the real noiseless case
    # mc_results = getBiasForBins(results_real,truth_array,bin_column=stats_array['rgp'], bin_ids=stats_array['id_cosmos'],bin_values=bins_redshift)
    # filename_results = 'bins.redshift.real.cat'
    # writeMCresults(filename_results,mc_results)



    # truth_array = tabletools.loadTable(table_name='truth_array',filepath=filepath_truth,dtype=dtype_table_truth,logger=logger)
    # results_array = tabletools.loadTable(table_name='results_array',filepath=filepath_results,dtype=dtype_table_results2,logger=logger)

    # # get total bias
    # logger.info('getting results for all')
    # analyse.getBiasForResults(results_array,truth_array,logger=logger,n_gals_per_mean=2e6)





def plotModelBias():

    def _binCenters(b):
          c1 = [ b[i-1] + (b[i] - b[i-1])/2. for i in range(1,len(b)) ] 
          return c1

    global results_stats

    file_pickle = open(args.filepath_stats)
    results_stats = pickle.load(file_pickle)
    logger.info('opened stats file %s with %d rows' % (args.filepath_stats,results_stats.shape[0]))


    filename_acs_join_stats = 'acs_joins_stats.pp'
    if os.path.isfile(filename_acs_join_stats):
        results_modd = tabletools.loadTable(table_name='results_modd',filepath=filename_acs_join_stats,logger=logger)
    else:
        results_modd = getTableACSjoinStats()
        # logger.info('getting modd')
        # acs_array = pyfits.getdata(args.filepath_acs,1)
        # results_modd = numpy.zeros(len(results_stats))
        # for i in range(len(results_modd)):
        #     select = numpy.nonzero(acs_array['IDENT']==results_stats['cosmos_id'][i])
        #     results_modd[i] = acs_array['modd'][select]
        # logger.info('got modd')
        # tabletools.saveTable(filename_acs_join_stats,results_modd)

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
    bins_snr = numpy.logspace(2,4,n_bins/2)
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
    bins_size = numpy.linspace(1,4,n_bins)
    h,b,p = pylab.hist(results_size,bins=bins_size)
    pylab.xlabel('rgp/rp')    
    filename_fig = 'figures/figure.hist.size.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # hlr distribution

    pylab.figure()
    results_hlr = results_stats['hlr']/pixel_scale
    bins_hlr = numpy.linspace(0,5,n_bins)
    h,b,p = pylab.hist(results_hlr,bins=bins_hlr)
    pylab.xlabel('hlr [pixels]')    
    filename_fig = 'figures/figure.hist.hlr.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # redshift distribution

    pylab.figure()
    bins_redshift = [0 , 0.35, 0.6, 0.8, 1.1 , 1.5]
    h,b,p = pylab.hist(results_zphot,bins=n_bins)
    # h,b,p = pylab.hist(results_zphot,bins=bins_redshift)
    pylab.xlabel('photo-z')    
    filename_fig = 'figures/figure.hist.redshift.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

    # modd distribution

    pylab.figure()
    # bins_modd = [s for s in set(results_modd)]
    bins_modd = [0,9,32]
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
    digitized = numpy.digitize(results_zphot, bins_redshift)
    for i in range(1, len(bins_redshift)):
        h,b = pylab.histogram(results_size[digitized == i],bins=bins_size,normed=True)
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
    colors_modd_m1 = ['r','b']
    colors_modd_m2 = ['m','c']
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

    import colorsys
    n_colors = len(bins_redshift)
    HSV_tuples = [(x*1.0/n_colors, 0.75, 0.75) for x in range(n_colors)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)

    for i in range(1, len(bins_redshift)):
        m1= results_m1[digitized == i]
        m2= results_m2[digitized == i]
        h,b,p = pylab.hist(m1,bins=n_bins_m,normed=True,histtype='step',label='m1,2 z=%2.2f' % (bins_redshift[i]),color=RGB_tuples[i-1])
        h,b,p = pylab.hist(m2,bins=n_bins_m,normed=True,histtype='step', color=RGB_tuples[i-1])
        pylab.axvline(linewidth=2, color='k')
        pylab.axvline(x=m1.mean(),linewidth=1, color=RGB_tuples[i-1])
        pylab.axvline(x=m2.mean(),linewidth=1, color=RGB_tuples[i-1])        

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
        h,b,p = pylab.hist(results_zphot[digitized == i],bins=bins_redshift,normed=True,histtype='step',label='type=%d' % i )
        
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
    for i in range(1, len(bins_modd)):
        # print set(digitized)
        h,b,p = pylab.hist(results_size[digitized == i],bins=bins_size,normed=True,histtype='step',label='type=%d' % i )
        
    pylab.grid()
    pylab.legend()
    pylab.xlabel('size')
    pylab.ylabel('normalised hist')

    filename_fig = 'figures/figure.hist.MMODforSIZE.real.png'
    pylab.savefig(filename_fig)
    logger.info('saved %s' % filename_fig)

def getBiasForBins(results_array,truth_array,bin_column, bin_ids,bin_values , bin_column_name='test'):

    """
    @brief plot bias for different bins of different parameters
    @results_array use this array to get the results -- it can be noiseless, noisy real or noisy bfit
    @truth_array truth array used to create results_array
    @info_string some info for plot titles
    """

    # initialise bins
    bins_ids = []
    digitized = numpy.digitize(bin_column, bin_values)
    bins_ids = [bin_ids[digitized == i] for i in range(1, len(bin_values))] 
    
    mc_results = []

    # loop over bins
    for i,ids in enumerate(bins_ids):
        results_bin,truth_bin,_ = analyse.selectByIDs(ids,results_array,truth_array,logger)
        logger.info('bin %d, number of galaxies in sample %5d' % (i,len(results_bin)))
        mcr = analyse.getBiasForResults(results_bin,truth_bin,logger=logger,bin_param=bin_column_name,bin_id=i)
        mcr['bin_id'] = i
        mcr['bin_value'] = bin_values[i]
        mc_results.append(mcr)

    mc_results = numpy.concatenate(mc_results)

    return mc_results

def writeMCresults(filename_results,mc_results):

    header = '# bin_id bin_value n_gals_in_bin m1 m2 c1 c2 std_m1 std_m2 std_c1 std_c2'
    fmt = '%d\t%f\t%d\t' + '%2.10e\t'*8 + '\n'

    file_results = open(filename_results,'w')
 
    for mcr in mc_results:
        line = fmt % (mcr['bin_id'],mcr['bin_value'],mcr['n_gals_in_bin'],mcr['m1'],mcr['m2'],mcr['c1'],mcr['c2'],mcr['std_m1'],mcr['std_m2'],mcr['std_c1'],mcr['std_c2'])
        file_results.write(line)

    file_results.close()
    logger.info('wrote file %s with %d bins' % (filename_results,len(mc_results)) )




    
def main():

    global logger , config , args

    description = 'Plots for the results of nmb_main. To use in ipython, create a variable global results_array, global truth_array to speed up loading'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
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

    # get total bias

    truth_array_25880  = tabletools.loadTable(table_name='truth_array_25880',  filepath = filepath_truth_25880,        dtype = dtype_table_truth,       logger=logger)
    truth_array_26000  = tabletools.loadTable(table_name='truth_array_26000',  filepath = filepath_truth_26000,        dtype = dtype_table_truth,       logger=logger)
    results_bfit_noisy = tabletools.loadTable(table_name='results_bfit_noisy', filepath = filepath_results_real_noisy, dtype = dtype_table_results2,    logger=logger)
    results_real_noisy = tabletools.loadTable(table_name='results_real_noisy', filepath = filepath_results_bfit_noisy, dtype = dtype_table_results2,    logger=logger)
    results_real       = tabletools.loadTable(table_name='results_real',       filepath = filepath_results_real,       dtype = dtype_table_results,     logger=logger)
    stats_array        = tabletools.loadTable(table_name='stats_array',        filepath = filepath_stats,              logger=logger)
    ajs_array          = tabletools.loadTable(table_name='ajs_array',          filepath = filepath_acs_join_stats,     logger=logger)
    acs_array          = tabletools.loadTable(table_name='acs_array',          filepath = filepath_acs,                logger=logger)

    mc_bfit_noisy = analyse.getBiasForResults(results_bfit_noisy,truth_array_25880,n_gals_per_mean=2e6,bin_param='all',bin_id=0)
    mc_real_noisy = analyse.getBiasForResults(results_real_noisy,truth_array_25880,n_gals_per_mean=2e6,bin_param='all',bin_id=0)

    print mc_bfit_noisy
    print mc_real_noisy

    # get the bins
    
    # saveBiasForBins()
    # plotBiasForBins()
    

if __name__ == "__main__":

    main()