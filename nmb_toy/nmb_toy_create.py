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
import i3_classes

sys.path.append("/home/kacprzak/code/tutils")

from mattools import *
from numpy import *
 
dirname_images = './images/'
dirname_ini = './ini/'
dirname_yaml = './yaml/'

filepath_galsim_yaml = '/home/kacprzak/code/GalSim/bin/galsim_yaml';
filename_ini_templ = 'nmb.ini.templ'
filename_yaml_templ = 'nmb_toy.yaml.templ'
filename_truth = 'truth.nmb_toy.txt'
filename_runlist = 'runlist.nmb_toy.txt'
filename_gals_cat = 'gals.cat'


def createRun():

	file_runlist = open(filename_runlist,'w')
	file_truth = open(filename_truth,'w')
	truth_header = '# true_si true_hlr true_g1 true_g2 model_si model_hlr model_g1 model_g2 noise_hlr noise_g1 noise_g2\n'
	file_truth.write(truth_header)
	
	grid_sersic_index = arange(config['grid']['min'],config['grid']['max']+config['grid']['step'],config['grid']['step'])

	file_conf_templ = open(filename_yaml_templ,'r')
	conf_templ = file_conf_templ.read()

	file_ini_templ = open(filename_ini_templ,'r')
	ini_templ = file_ini_templ.read()

	logging.info('got %d galaxies' % len(config['gals']))

	for iser,ser in enumerate(grid_sersic_index):

		ini_filled = ini_templ % (config['n_pix'],ser)
		filename_ini = 'si%02d.ini' % iser
		filepath_ini = os.path.join(dirname_ini,filename_ini)
		file_ini = open(filepath_ini,'w')
		file_ini.write(ini_filled)
		file_ini.close()

	for igal,gal in enumerate(config['gals']):

		# get the real galaxy
		n_tile=config['n_tile']
		n_pix=config['n_pix']
		hlr = gal['half_light_radius']
		g1  = gal['g1']
		g2  = gal['g2']
		sersic_index_real = gal['sersic_index']
		filename_real = 'real%02d.fits' % (igal)
		filepath_real = os.path.join(dirname_images,filename_real)

		conf_filled = conf_templ % (sersic_index_real,hlr,g1,g2,n_tile,n_tile,n_pix,n_pix,filepath_real)

		filename_conf = 'real%02d.yaml' % igal
		filepath_conf = os.path.join(dirname_yaml,filename_conf)
		file_conf = open(filepath_conf,'w')
		file_conf.write(conf_filled)
		file_conf.close()

		filename_cmd = 'cmd.sh'
		file_cmd = open(filename_cmd,'w')
		file_cmd.write('python %s %s\n' % (filepath_galsim_yaml, filepath_conf))
		file_cmd.close()

		if (not args.reimage) and os.path.isfile(filepath_real):
			logging.info('NOT creating %s' % filepath_real)
		else:
			subprocess.call(('sh',filename_cmd))

		# filename_ini = os.path.join(dirname_ini,'si%02d.ini' % ser)

		image_tiled = pyfits.getdata(filepath_real)
		image_stamp = image_tiled[0:n_pix,0:n_pix]
		noise_std = linalg.norm(image_stamp)/config['snr']
		img_real = image_tiled
		noise_same = []

		for nn in range(config['n_reps_diff']):

			noise = random.randn(img_real.shape[0],img_real.shape[1])*noise_std
			filename_noisy_real_diff = filename_real + ('.d%02d' % nn)
			filepath_noisy_real_diff = os.path.join(dirname_images,filename_noisy_real_diff)
			img_real_noisy_diff = img_real + noise

			if (not args.reimage) and os.path.isfile(filepath_noisy_real_diff):
				logging.info('NOT creating noisy images %s' % (filepath_noisy_real_diff))
			else:
				pyfits.writeto(filepath_noisy_real_diff,img_real_noisy_diff.astype(float32),clobber=True)

			for iser,ser in enumerate(grid_sersic_index):	

				filename_ini = 'si%02d.ini' % iser
				filepath_ini = os.path.join(dirname_ini,filename_ini)

				file_runlist.write('%s\t%s\n' % (filename_noisy_real_diff,filename_ini))

		for nn in range(config['n_reps_same']):

			noise_same.append(random.randn(img_real.shape[0],img_real.shape[1])*noise_std)
			filename_noisy_real_same = filename_real + ('.s%02d' % nn)
			filepath_noisy_real_same = os.path.join(dirname_images,filename_noisy_real_same)
			img_real_noisy_same = img_real + noise_same[nn]

			if (not args.reimage) and os.path.isfile(filepath_noisy_real_same):
				logging.info('NOT creating noisy images %s' % (filepath_noisy_real_same))
			else:
				pyfits.writeto(filepath_noisy_real_same,img_real_noisy_same.astype(float32),clobber=True)

			for iser,ser in enumerate(grid_sersic_index):	

				filename_ini = 'si%02d.ini' % iser
				filepath_ini = os.path.join(dirname_ini,filename_ini)

				file_runlist.write('%s\t%s\n' % (filename_noisy_real_same,filename_ini))



		# now create the bestfit images
		for iser,ser in enumerate(grid_sersic_index):

			# this should have been created already
			filename_ini = 'si%02d.ini' % iser
			filepath_ini = os.path.join(dirname_ini,filename_ini)
			logging.info('running im3shape')
			i3gal = getBestFit(image_stamp,filepath_ini)
			ser_g1 = i3gal.params_gal_measured[2]
			ser_g2 = 0
			ser_hlr  = i3gal.params_gal_measured[4]

			truth_line = '%2.2f\t%2.8f\t% 2.8f\t% 2.8f\t%2.2f\t%2.8f\t% 2.8f\t% 2.8f\n' % (sersic_index_real,hlr,g1,g2,ser,ser_hlr,ser_g1,ser_g2)
			logging.info(truth_line)
			file_truth.write(truth_line)

			filename_bfit = 'real%02d.bfit%02d.fits' % (igal,iser)
			filepath_bfit = os.path.join(dirname_images,filename_bfit)

			conf_filled = conf_templ % (ser,ser_hlr,ser_g1,ser_g2,n_tile,n_tile,n_pix,n_pix,filepath_bfit)

			filename_conf = 'real%02d.bfit%02d.yaml' % (igal,iser)
			filepath_conf = os.path.join(dirname_yaml, filename_conf)
			file_conf = open(filepath_conf,'w')
			file_conf.write(conf_filled)
			file_conf.close()

			filename_cmd = 'cmd.sh'
			file_cmd = open(filename_cmd,'w')
			file_cmd.write('python %s %s\n' % (filepath_galsim_yaml, filepath_conf))
			file_cmd.close()

			if (not args.reimage) and os.path.isfile(filepath_bfit):
				logging.info('NOT creating %s' % filepath_bfit)
			else:
				subprocess.call(('sh',filename_cmd))

			# create nosiy versions

			noise_std = linalg.norm(image_stamp)/config['snr']
	
			img_real = pyfits.getdata(filepath_real)
			img_bfit = pyfits.getdata(filepath_bfit)

			for nn in range(config['n_reps_diff']):

				# add different noise maps
				filename_noisy_bfit_diff = filename_bfit + ('.d%02d' % nn)
				filepath_noisy_bfit_diff = os.path.join(dirname_images,filename_noisy_bfit_diff)
				noise = random.randn(img_bfit.shape[0],img_bfit.shape[1])*noise_std
				img_bfit_noisy_diff = img_bfit + noise
				
				if (not args.reimage) and os.path.isfile(filepath_noisy_bfit_diff):
					logging.info('NOT creating noisy images %s' % (filepath_noisy_bfit_diff))
				else:
					pyfits.writeto(filepath_noisy_bfit_diff,img_bfit_noisy_diff.astype(float32),clobber=True)

				file_runlist.write('%s\t%s\n' % (filename_noisy_bfit_diff,filename_ini))

			for nn in range(config['n_reps_same']):
					
				# add same noise maps
				filename_noisy_bfit_same = filename_bfit + ('.s%02d' % nn)
				filepath_noisy_bfit_same = os.path.join(dirname_images,filename_noisy_bfit_same)
				img_bfit_noisy_same = img_bfit + noise_same[nn]

				if (not args.reimage) and os.path.isfile(filepath_noisy_bfit_same):
					logging.info('NOT creating noisy image %s' % (filepath_noisy_bfit_same))
				else:
					pyfits.writeto(filepath_noisy_bfit_same,img_bfit_noisy_same.astype(float32),clobber=True)

				file_runlist.write('%s\t%s\n' % (filename_noisy_bfit_same,filename_ini))

def getBestFit(image_obj,file_ini="nmb.create.ini"):
	
	i3gal = i3_classes.Galaxy()
	i3gal.load_ini_file(file_ini)
	i3gal.params_gal_true[0] = image_obj.shape[0]/2.
	i3gal.params_gal_true[1] = image_obj.shape[1]/2.
	
	# print image_obj.array.shape

	i3gal.image_noisy = image_obj.copy()
	# i3gal.image_psf_hires = image_psf_hires.copy()
	# change this if sersic_templ changes
	i3gal.params_psf_true = array([3., 2.85 ,0., 0.])
			
	i3gal.load_ini_file(file_ini=file_ini)
	i3gal.get_measurement(show_cmd=False)
	
	logging.debug(arr2str(i3gal.params_gal_measured))
	
	return i3gal


def getGalsCatalog():
	
	cat_gal = open(filename_gals_cat,'w')

	id = 0

	for ix1 in range(config['n_tile']):
		
		for ix2 in range(config['n_tile']):
			
			id = id + 1
			
			bx1 = ix1 * config['n_pix'] + config['n_pix']/2.
			bx2 = ix2 * config['n_pix'] + config['n_pix']/2.
			
			row = "%d %2.1f %2.1f \n" % (id,bx1,bx2)
			
			cat_gal.write(row)
			# cat_psf.write(psf_info)
			
	logging.info("saved catalogue files %s" % filename_gals_cat)
		

description = """
in progress
"""

parser = argparse.ArgumentParser(description=description, add_help=True)
parser.add_argument('filename_config', type=str, help='yaml config file, see example_config.yaml for example.')
parser.add_argument('--reimage', action="store_true", help='Whether get the images again even if they exist', default=False)


global args
args = parser.parse_args()

filename_config = args.filename_config
global config
config = yaml.load(open(filename_config,'r'))

# global real_galaxy_catalog

format_string= "%(asctime)-15s %(message)s"
logging.basicConfig(filename=args.filename_config + '.log',level=logging.INFO,format=format_string)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

createRun()
getGalsCatalog()