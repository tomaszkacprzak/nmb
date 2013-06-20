import os
import pdb
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
import tabletools
from tablespec import *

def getFWHM(i3_result,fwxm=0.5,n_sub=3):


        n_pix = config['image']['size']*n_sub
        pixel_scale           = config['image']['pixel_scale']
        pixel_scale_upsampled = config['image']['pixel_scale']/float(n_sub)
        gal1 = galsim.DeVaucouleurs(half_light_radius=i3_result.sersic_parameter_radius*pixel_scale)
        gal2 = galsim.Exponential(half_light_radius=i3_result.sersic_parameter_radius*pixel_scale)
        gal = i3_result.sersic_disc_flux*gal2 + i3_result.sersic_bulge_flux*gal1
        psf_beta = config['psf']['beta']
        psf_fwhm = config['psf']['fwhm']
        psf = galsim.Moffat(beta=psf_beta,fwhm=psf_fwhm)
        pix = galsim.Pixel(xw=pixel_scale)
        obj = galsim.Convolve([gal,psf,pix])
        img = galsim.ImageD(n_pix,n_pix)
        obj.draw(img,dx=pixel_scale_upsampled)

        # fhwm code

        image_nx = n_pix
        # xy0 = n_pix/2 + 0.5 + shift;
        xy0 = n_pix/2 + 0.5;

        # set the minimum element to 0
        image_true = img.array
        image_true = image_true - numpy.amin(image_true.flatten()) 
                
        f1 = 0.
        f2 = 0.
        x1 = 0
        x2 = 0 
                
        profile = image_true[int(image_nx//2),:]
                
        max_ind = int(image_nx//2)
        max_val = profile[max_ind]
        f3 = max_val*fwxm
        
        diff = abs(profile-f3)
        
        x1 = numpy.argmin(diff)
        f1 = profile[x1]
        
        if( f1 < f3 ):  x2 = x1+1
        else:       x2 = x1-1
        f2 = profile[x2];
    
        a = (f1-f2)/(x1 - x2)
        b = f1 - a*x1; 
        x3 = (f3 - b)/a;
               
        fwhm = 2.*abs(max_ind-x3) * pixel_scale_upsampled                         
        return fwhm

def getGalaxyImages():

    config1 = copy.deepcopy(config)

    # get the total number of objects
    n_obj = len(galsim.config.GetNObjForMultiFits(config1,0,0))

    logger.info('found %d objects in the config file' % n_obj)

    # make sure that we are using correct obj_num and nimages
    if config1['args'].obj_num < 0: raise ValueError('start_id should be non-negative')
    else: obj_num = config1['args'].obj_num   
    if config1['args'].nimages <= 0: nimages = n_obj
    else: nimages = config1['args'].nimages

    logger.info('starting with obj_num=%6d and processing %6d objects' % (obj_num,nimages))    
    
    # process input  
    galsim.config.ProcessInput(config1)

    logger.info('getting images')

    # get the galsim galaxy
    img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config1,obj_num=obj_num,nimages=nimages,make_psf_image=True,logger=logger_config)
    logger.info('got %d images' % len(img_gal))

    return img_gal

def runIm3shape():

    # open the RGC
    if 'real_catalog' in config['input']:
        filepath_rgc = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
        rgc = pyfits.open(filepath_rgc)[1]

    # open the ring test catalog
    filename_cat = os.path.join(config['input']['catalog']['dir'],config['input']['catalog']['file_name'])
    # truth_cat = numpy.loadtxt(filename_cat,dtype=dtype_table_truth)
    truth_cat = tabletools.loadTable(filepath=filename_cat,table_name='truth_cat',dtype=dtype_table_truth,logger=logger)

    n_objects = truth_cat.shape[0]

    # get im3shape
    dirpath_im3shape = os.path.join(os.environ['IM3SHAPE'],'python')
    sys.path.append(dirpath_im3shape)
    import im3shape

    # get images
    img_gals = getGalaxyImages()

    # get options
    n_pix = config['image']['size']
    pixel_scale = config['image']['pixel_scale']
    i3_options = im3shape.I3_options()
    i3_options.read_ini_file(config['args'].filepath_ini)
    logger.info('loaded im3shape ini file %s' % config['args'].filepath_ini)  
    i3_options.stamp_size = n_pix

    # get the file 
    filename_results = 'results.%s.%012d.cat' % (config['args'].name_config,config['args'].obj_num)
    file_results = open(filename_results,'w')

    # create PSF from Moffat parameters
    # get i3 images - get first i3_galaxy to initialise the PSF - kind of strange, but hey..
    i3_galaxy = im3shape.I3_image(n_pix, n_pix)
    psf_beta = float(config['psf']['beta'])
    psf_fwhm = float(config['psf']['fwhm'])/float(pixel_scale)
    psf_e1 = float(config['psf']['ellip']['g1'])
    psf_e2 = float(config['psf']['ellip']['g2'])
    i3_psf = i3_galaxy.make_great10_psf(psf_beta, psf_fwhm, psf_e1, psf_e2, i3_options)

    obj_num = config['args'].obj_num

    # loop over all created images
    for ig,img_gal in enumerate(img_gals):

        # get i3 images
        i3_galaxy = im3shape.I3_image(n_pix, n_pix)
        scale = img_gal.array.sum()
        i3_galaxy.from_array(img_gal.array)  
        i3_galaxy.scale(1.0/scale)


        # this is a workaround - the array for model bias real images had id_unique,
        # but the results files have 'identifier' - so we use one or the other id they exist
        # get the unique_id
        if 'id_unique' in truth_cat.dtype.names:

            id_global = ig
            id_object = (obj_num+ig) % n_objects
            id_cosmos = truth_cat['id_cosmos'][ id_object ]
            id_unique = truth_cat['id_unique'][ id_object ]

        elif 'identifier' in truth_cat.dtype.names:
            
            id_global = ig
            id_object = (obj_num+ig) % n_objects
            id_unique = truth_cat['identifier'][ id_object ]
            id_cosmos = int(truth_cat['identifier'][ id_object ]//10000)

        
        i3_result, i3_best_fit = im3shape.i3_analyze(i3_galaxy, i3_psf, i3_options, ID=id_unique)


        saveResult(file_results,i3_result,id_global,id_object,id_unique,id_cosmos)
        printResult(i3_result,id_global)
        if 'e1' in truth_cat.dtype.names: printTruth(i3_result,truth_cat[ig])


        # save residual plots
        if config['args'].verbosity > 2:

            i1 = i3_best_fit.array/sum(i3_best_fit.array.flatten())
            i2 = img_gal.array/sum(img_gal.array.flatten())          

            import pylab
            pylab.subplot(1,5,1)
            pylab.imshow(i1,interpolation='nearest')
            pylab.title('best fit')

            pylab.subplot(1,5,2)
            pylab.imshow(i2,interpolation='nearest')
            pylab.title('galaxy')

            pylab.subplot(1,5,3)
            pylab.imshow(i1-i2,interpolation='nearest')
            pylab.title('residuals')

            pylab.subplot(1,5,4)
            pylab.imshow(i3_psf.array,interpolation='nearest')
            pylab.title('PSF')

            pylab.subplot(1,5,5)
            pylab.imshow(img_gal.array,interpolation='nearest')
            pylab.colorbar()
            pylab.title('best fit')


            filename_fig = 'debug/fig.residual.%09d.png' % id_unique
            pylab.savefig(filename_fig)
            pylab.close()


def saveResult(file_results,i3_result,idg,ido,idu,idc):

    pixel_scale = config['image']['pixel_scale']
    n_pix = config['image']['size']

    fmt = '%d\t%d\t%d\t%d\t% e\t% 2.2f\t' + '% e\t'*5 + '%2.2f\t' + '% e\t'*8 + '% 5d'*3 + '\n'
    line = fmt % (
                 idg,
                 ido,
                 idu,
                 idc,
                 i3_result.likelihood,
                 i3_result.time_taken,
                 (i3_result.sersic_parameter_x0 - float(n_pix)/2.)*pixel_scale*(-1.), # wierd flip, see compare_im3shape_galsim
                 (i3_result.sersic_parameter_y0 - float(n_pix)/2.)*pixel_scale*(-1.),
                 i3_result.sersic_parameter_e1,
                 i3_result.sersic_parameter_e2,
                 i3_result.sersic_parameter_radius*pixel_scale,
                 getFWHM(i3_result),
                 i3_result.sersic_bulge_flux,
                 i3_result.sersic_disc_flux,
                 i3_result.sersic_flux_ratio,
                 i3_result.stats_signal_to_noise,
                 i3_result.stats_min_residuals,
                 i3_result.stats_max_residuals,
                 i3_result.stats_model_min,
                 i3_result.stats_model_max,
                 i3_result.levmar_number_of_likelihood_evaluations,
                 i3_result.levmar_number_of_iterations,
                 i3_result.levmar_reason_of_termination,
                 )
    
    file_results.write(line)

def printResult(i3_result,id_global):


    pixel_scale = config['image']['pixel_scale']
    n_pix = config['image']['size']

    fmt = '%10d\t%10d\t% e\t% 2.2f\t' + '% e\t'*5 + '%2.2f\t'*2
    line = fmt % (
                 id_global,
                 i3_result.identifier,
                 i3_result.likelihood,
                 i3_result.time_taken,
                 i3_result.sersic_parameter_e1,
                 i3_result.sersic_parameter_e2,
                 i3_result.sersic_parameter_radius*pixel_scale,
                 i3_result.sersic_bulge_flux,
                 i3_result.sersic_disc_flux,
                 i3_result.stats_signal_to_noise,
                 getFWHM(i3_result),
                 )
    
    logger.info(line)

def printTruth(i3_result,truth_row):

    pixel_scale = config['image']['pixel_scale']
    n_pix = config['image']['size']

    fmt = '%d\t% e\t% 2.2f\t' + '% e\t'*5
    line = fmt % (
                 i3_result.identifier,
                 truth_row['likelihood'],
                 truth_row['time_taken'],
                 truth_row['e1'],   
                 truth_row['e2'],  
                 truth_row['radius']*pixel_scale ,
                 truth_row['bulge_flux'],
                 truth_row['disc_flux']
                 )
    
    logger.info(line)

if __name__ == "__main__":

    description = 'Noise and model bias driver. Requires $IM3SHAPE and $GALSIM to be set.'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('filepath_config', type=str, help='yaml config file, see nmb_main.real.test.yaml for example.')
    parser.add_argument('--filepath_ini', type=str, default='nmb.ini', help='ini im3shape config file (default: nmb.ini)')
    parser.add_argument('--filepath_truth', type=str, default=None, help='truth file for the run, overrides the config file (by default is taken from yaml file)')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-snr', '--signal_to_noise', type=float, action='store', default=None, help='signal to noise at which to run the test')
    parser.add_argument('--obj_num',  type=int, action='store', default= 0, help= 'first obj_num in config to process (starts from 1)') 
    parser.add_argument('--nimages',  type=int, action='store', default=-1, help= 'number of images to process, starting with obj_num')
    
    args = parser.parse_args()
    args.name_config = os.path.basename(args.filepath_config).replace('.yaml','')

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    global logger , logger_config , config
    logger = logging.getLogger("nmb_main") 
    config_logging_level = logging_levels[args.verbosity - 1]
    logger_config = logging.getLogger("galsim_config") 
    logger_config.setLevel(config_logging_level)

    if args.verbosity == 3: 
        if not os.path.exists('./debug'): 
            os.makedirs('./debug')

    # load the configuration file
    config = yaml.load(open(args.filepath_config,'r')) 
    # store the args in config so it's easier to use them
    config['args'] = args
    # change config to match signal to noise
    if args.signal_to_noise != None: config['gal']['signal_to_noise'] = args.signal_to_noise
           
    # load site config
    if 'real_catalog' in config['input']:
        config['input']['real_catalog']['dir'] = os.path.join(os.environ['GALSIM'],'rgc')
        logger.info('rgc set to be in %s' % config['input']['real_catalog']['dir'])

    # set the truth file path
    if args.filepath_truth != None:
            config['input']['catalog']['dir'] = os.path.dirname(args.filepath_truth)
            config['input']['catalog']['file_name'] = os.path.basename(args.filepath_truth)
            logger.info('truth catalog set to be in %s %s' % (config['input']['catalog']['dir'],config['input']['catalog']['file_name']))

    # run the main driver
    runIm3shape()



