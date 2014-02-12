import tabletools
import pylab as pl
import scipy.stats
import numpy as np

import plotstools




cosmos=tabletools.loadTable('cosmos_acs_shera_may2011.fits.gz')
cos=cosmos[(cosmos['FWHM_IMAGE']<100) * (cosmos['ZPHOT']>0) * (cosmos['ZPHOT']<3) ]

npix = 50
pixel_size=0.05
plotstools.plot_dist(X=[cos['FWHM_IMAGE']*pixel_size,cos['ZPHOT'],cos['MODD']],bins=[30,30,range(0,31)],labels=['FWHM [arcsec]' , 'ZPHOT' , 'MODD'])
filename_fig = 'histograms_size_zphot_modd.eps'

left  = 0.1  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.0   # the amount of width reserved for blank space between subplots
hspace = 0.0   # the amount of height reserved for white space between subplots
pl.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
pl.savefig(filename_fig)


print 'saved' , filename_fig





#like is 2D, the axes are 1D
