import glob
import shutil
import yaml
import random
import os
import pyfits
import argparse
import subprocess
import sys
import galsim
import copy
import logging
import pylab
import scipy
import scipy.io
import pdb
import i3_classes

sys.path.append("/home/kacprzak/code/tutils")

from mattools import *
from numpy import *
 

dirname_images = './images/'
dirname_results = './results/'
dirname_ini = './ini/'
dirname_yaml = './yaml/'
dirname_results_mat = './results_mat/'

filepath_galsim_yaml = '/home/kacprzak/code/GalSim/bin/galsim_yaml';
filename_ini_templ = 'nmb.ini.templ'
filename_yaml_templ = 'nmb_toy.yaml.templ'
filename_truth = 'nmb_toy.truth.txt'
filename_stats = 'nmb_toy.stats.%s.txt'
filename_runlist = 'runlist.nmb_toy.txt'
filename_result_tmp = 'tmp.result.txt'

n_same_reals_use = 5e3

def mergeReal(filename_img,filename_ini,noise_type):


	filename_result_mat = os.path.join(dirname_results_mat,'%s.%s.%s.nbc.mat' % (filename_img,noise_type,filename_ini)) 

	if not os.path.isfile(filename_result_mat):

		filename_result_templ = os.path.join(dirname_results,'%s.%s*.%s.nbc*' % (filename_img,noise_type,filename_ini)) 
		files_found = glob.glob(filename_result_templ)
		n_parts = len(files_found)
		results = []
		for i,f in enumerate(files_found):
			# ./results/real00.fits.d29.si00.ini.nbc.436  36
			part_id = int(f[23:25])
			# print f + '  ' + str(part_id) 
			results_part = loadtxt(f)
			# add an id*1000 to be able to match up noises in case there are some parts missing
			results_part[:,0] = results_part[:,0] + part_id*1e4

			if i==0:
				results = results_part
			else:
				results = concatenate([results,results_part]);


		matdict = {}
		matdict['results']=results
		
		scipy.io.savemat(filename_result_mat,matdict)
		
		logging.info('getting galaxy from image %s - found %d parts  - created file %s' % (filename_result_templ,n_parts,filename_result_mat))

	else:
		
		logging.info('getting galaxy from image %s - file %s already exists ' % (filename_img,filename_result_mat))

 
 
def mergeBfit(filename_img,filename_ini,noise_type):


	filename_result_mat = os.path.join(dirname_results_mat,'%s.%s.%s.nbc.mat' % (filename_img,noise_type,filename_ini)) 

	if not os.path.isfile(filename_result_mat):

		filename_result_templ = os.path.join(dirname_results,'%s.%s*.%s.nbc.*' % (filename_img,noise_type,filename_ini)) 
		files_found = glob.glob(filename_result_templ)
		n_parts = len(files_found)
		results = []
		for i,f in enumerate(files_found):
			part_id = int(f[30:32])
			# print f + '  ' + str(part_id) 
			results_part = loadtxt(f)
			# add an id*1000 to be able to match up noises in case there are some parts missing
			results_part[:,0] = results_part[:,0] + part_id*1e4

			if i==0:
				results = results_part
			else:
				results = concatenate([results,results_part]);

		matdict = {}
		matdict['results']=results
		
		scipy.io.savemat(filename_result_mat,matdict)
		


		logging.info('getting galaxy from image %s - found %d parts  - created file %s' % (filename_result_templ,n_parts,filename_result_mat))

	else:
		
		logging.info('getting galaxy from image %s - file %s already exists ' % (filename_img,filename_result_mat))




def getStats(noise_type):

	logging.info('got %d galaxies' % len(config['gals']))

	grid_sersic_index = arange(config['grid']['min'],config['grid']['max']+config['grid']['step'],config['grid']['step'])

	truth_header = '# true_si true_hlr true_g1 true_g2 model_si model_hlr model_g1 model_g2 noise_hlr noise_g1 noise_g2 noise_g1_std noise_g2_std\n'

	stats_truth = loadtxt(filename_truth)
	n_points = stats_truth.shape[0]
	n_entries = stats_truth.shape[0]

	stats_resutls = concatenate([stats_truth,zeros([n_points,20])],axis=1)

	irow = 0

	for igal,gal in enumerate(config['gals']):
				for iser,ser in enumerate(grid_sersic_index):

# get stats for bfit

					filename_result_bfit_mat = os.path.join(dirname_results_mat,'real%02d.bfit%02d.fits.%s.si%02d.ini.nbc.mat' % (igal,iser,noise_type,iser)) 			
					filename_result_real_mat = os.path.join(dirname_results_mat,'real%02d.fits.%s.si%02d.ini.nbc.mat' % (igal,noise_type,iser)) 

					mat = scipy.io.loadmat(filename_result_bfit_mat); results_bfit = mat['results']

					if results_bfit.size == 0: results_bfit=ones([666,24])						

					n_gals = results_bfit.shape[0]
		
					results_mean = mean(results_bfit,axis=0)
					results_stdv = std(results_bfit,axis=0,ddof=1)

					mean_e1 = results_mean[cols['results']['e1']]
					mean_e2 = results_mean[cols['results']['e2']]
					mean_re = results_mean[cols['results']['re']]

					stdv_e1 = results_stdv[cols['results']['e1']]
					stdv_e2 = results_stdv[cols['results']['e2']]
					stdv_re = results_stdv[cols['results']['re']]

					stdm_e1 = stdv_e1/sqrt(float(n_gals)) 
					stdm_e2 = stdv_e2/sqrt(float(n_gals)) 
					stdm_re = stdv_re/sqrt(float(n_gals)) 

					stats_resutls[irow,cols['stats']['noise_bfit_hlr']] = mean_re
					stats_resutls[irow,cols['stats']['noise_bfit_g1']]  = mean_e1
					stats_resutls[irow,cols['stats']['noise_bfit_g2']]  = mean_e2
					stats_resutls[irow,cols['stats']['noise_bfit_hlr_stdm']] = stdm_re
					stats_resutls[irow,cols['stats']['noise_bfit_g1_stdm']]  = stdm_e1
					stats_resutls[irow,cols['stats']['noise_bfit_g2_stdm']]  = stdm_e2

					logging.info('%-50s\tbfit n_gals=%6d\t\tmean=[% 2.4f % 2.4f]\tstdv=[% 2.4f % 2.4f]\tstdm=[% 2.4f % 2.4f]' % (filename_result_bfit_mat, results_bfit.shape[0], mean_e1, mean_e2, stdv_e1, stdv_e2, stdm_e1, stdm_e2))					

# get stats for real

					mat = scipy.io.loadmat(filename_result_real_mat); results_real = mat['results']

					if results_real.size == 0: results_real=ones([666,24])						

					n_gals = results_real.shape[0]
		
					results_mean = mean(results_real,axis=0)
					results_stdv = std(results_real,axis=0,ddof=1)

					mean_e1 = results_mean[cols['results']['e1']]
					mean_e2 = results_mean[cols['results']['e2']]
					mean_re = results_mean[cols['results']['re']]

					stdv_e1 = results_stdv[cols['results']['e1']]
					stdv_e2 = results_stdv[cols['results']['e2']]
					stdv_re = results_stdv[cols['results']['re']]

					stdm_e1 = stdv_e1/sqrt(float(n_gals)) 
					stdm_e2 = stdv_e2/sqrt(float(n_gals)) 
					stdm_re = stdv_re/sqrt(float(n_gals)) 

					stats_resutls[irow,cols['stats']['noise_real_hlr']] = mean_re
					stats_resutls[irow,cols['stats']['noise_real_g1']]  = mean_e1
					stats_resutls[irow,cols['stats']['noise_real_g2']]  = mean_e2
					stats_resutls[irow,cols['stats']['noise_real_hlr_stdm']] = stdm_re
					stats_resutls[irow,cols['stats']['noise_real_g1_stdm']]  = stdm_e1
					stats_resutls[irow,cols['stats']['noise_real_g2_stdm']]  = stdm_e2

					logging.info('%-50s\treal n_gals=%6d\t\tmean=[% 2.4f % 2.4f]\tstdv=[% 2.4f % 2.4f]\tstdm=[% 2.4f % 2.4f]' % (filename_result_real_mat, results_real.shape[0], mean_e1, mean_e2, stdv_e1, stdv_e2, stdm_e1, stdm_e2))					


					# get the diff stats
	
					if noise_type == 's':

						n_use = n_same_reals_use

						results_bfit_common = results_bfit[in1d(results_bfit[:,0],results_real[:,0]),:]
						results_real_common = results_real[in1d(results_real[:,0],results_bfit[:,0]),:]

						results_bfit_common = results_bfit_common[argsort(results_bfit_common[:,0]),:]
						results_real_common = results_real_common[argsort(results_real_common[:,0]),:]

						results_bfit_common = results_bfit_common[0:n_use,:]
						results_real_common = results_real_common[0:n_use,:]

						if results_bfit_common.shape != results_real_common.shape:
							results_diff = ones([666,24])
						else:
							results_diff = results_real_common - results_bfit_common

						# pdb.set_trace()

						results_mean = mean(results_diff,axis=0)
						results_stdv = std(results_diff,axis=0,ddof=1)

						mean_e1 = results_mean[cols['results']['e1']]
						mean_e2 = results_mean[cols['results']['e2']]
						mean_re = results_mean[cols['results']['re']]

						stdv_e1 = results_stdv[cols['results']['e1']]
						stdv_e2 = results_stdv[cols['results']['e2']]
						stdv_re = results_stdv[cols['results']['re']]

						stdm_e1 = stdv_e1/sqrt(float(n_use)) 
						stdm_e2 = stdv_e2/sqrt(float(n_use)) 
						stdm_re = stdv_re/sqrt(float(n_use)) 

						stats_resutls[irow,cols['stats']['diff_hlr']] = mean_re
						stats_resutls[irow,cols['stats']['diff_g1']]  = mean_e1
						stats_resutls[irow,cols['stats']['diff_g2']]  = mean_e2
						stats_resutls[irow,cols['stats']['diff_hlr_stdm']] = stdm_re
						stats_resutls[irow,cols['stats']['diff_g1_stdm']]  = stdm_e1
						stats_resutls[irow,cols['stats']['diff_g2_stdm']]  = stdm_e2
						stats_resutls[irow,cols['stats']['diff_g1_stdv']]  = stdv_e1
						stats_resutls[irow,cols['stats']['diff_g2_stdv']]  = stdv_e2


						logging.info('%-50s\tdiff n_gals=%6d\t\tmean=[% 2.4f % 2.4f]\tstdv=[% 2.4f % 2.4f]\tstdm=[% 2.4f % 2.4f]' % (filename_result_real_mat, results_diff.shape[0], mean_e1, mean_e2, stdv_e1, stdv_e2, stdm_e1, stdm_e2))					

						n_bins=100
						bins = linspace(-0.2,0.2,n_bins)
						pylab.clf()
						pylab.hist(results_diff[:,4],bins)
						filename_fig = 'figures/hist.diff.real%02d.bfit%02d.png' % (igal,iser)
						pylab.savefig(filename_fig,format='png')


					irow = irow + 1

# get the difference

	filename = filename_stats % noise_type
	savetxt(filename,stats_resutls,fmt='% 2.4e')
	logging.info('wrote stats of results to %s',filename)




def mergeAll(noise_type):

	logging.info('got %d galaxies' % len(config['gals']))

	grid_sersic_index = arange(config['grid']['min'],config['grid']['max']+config['grid']['step'],config['grid']['step'])

	# merge real fitted with different sersic indices
	for igal,gal in enumerate(config['gals']):
		for iser,ser in enumerate(grid_sersic_index):

			# merge the real
			filename_img = 'real%02d.fits' % igal
			filename_ini = 'si%02d.ini' % iser
			mergeReal(filename_img,filename_ini,noise_type)

	# merge bfit fitted with same sersic index
	for igal,gal in enumerate(config['gals']):
		for iser,ser in enumerate(grid_sersic_index):

			# merge the bfit
			filename_img = 'real%02d.bfit%02d.fits' % (igal,iser)
			filename_ini = 'si%02d.ini' % iser
			mergeBfit(filename_img,filename_ini,noise_type)

description = """
in progress
"""

parser = argparse.ArgumentParser(description=description, add_help=True)
parser.add_argument('command', type=str, help='what to do? use {merge,stats}')
parser.add_argument('filename_config', type=str, help='yaml config file, see example_config.yaml for example.')
parser.add_argument('--filename_columns', action='store', type=str, help='yaml file, containing the column information of tables used.',default='tables_info.yaml')


global args
args = parser.parse_args()

filename_config = args.filename_config
global config
config = yaml.load(open(filename_config,'r'))


global cols
cols = yaml.load(open(args.filename_columns,'r'))


# global real_galaxy_catalog

format_string= "%(asctime)-15s %(message)s"
logging.basicConfig(filename=args.filename_config + '.analyse.log',level=logging.INFO,format=format_string)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

# analyseRun()
# getGalsCatalog()

if args.command == 'merge':
	if not os.path.exists('results_mat'): os.makedirs('results_mat')
	logging.info('----------------------------------------')
	logging.info('merging for different noise maps')
	logging.info('----------------------------------------')
	mergeAll('d')
	logging.info('----------------------------------------')
	logging.info('merging for same noise maps')
	logging.info('----------------------------------------')
	mergeAll('s')
elif args.command == 'stats':
	logging.info('----------------------------------------')
	logging.info('stats for different noise maps')
	logging.info('----------------------------------------')
	getStats('d')
	logging.info('----------------------------------------')
	logging.info('stats for same noise maps')
	logging.info('----------------------------------------')
	getStats('s')
else:
	logger.error('unknown command')
	raise ValueError('unknown command')
