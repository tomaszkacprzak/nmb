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

def getRingID(id_angle,id_shear):
    return id_shear*100 + id_angle

def getUniqueID(id_object,id_angle,id_shear):
    return id_object*10000 + getRingID(id_angle,id_shear)

def getUniqueID2(id_object,id_ring):
    return id_object*10000 + id_ring

def saveRingTestCatalog():

    def _writeShears(n_angles,shears,filename_cat):

        file_cat = open(filename_cat,'w')
        file_cat.write('# g1 g2 rotation_angle\n')
        fmt = '%04d\t% 2.4e\t% 2.4e\t% 2.10e\t%2d\t%2d\n' 
        d_angle = 180./n_angles    
        for ig,g in enumerate(shears):
            shears_g12 = [ [ 0.0, +g ] , [ 0.0, -g ] , [ +g,   0.0 ] , [ -g,   0.0 ] , [ +g,  +g ] , [-g,  -g ] , [ +g,  -g ] , [ -g,  +g ]]
            for ig12,g12 in enumerate(shears_g12):
                for ia in range(n_angles):
                    angle = ia*d_angle
                    id_ring = getRingID(ia,8*ig+ig12) 
                    line = fmt % (id_ring, g12[0], g12[1],   angle, ia, ig12); file_cat.write(line)
        file_cat.close()

    # this is for the main run
    # set number of angles in ring test
    n_angles = 8
    # set shears magnitutes, 8 shear for each mag will be created
    shears = [0.05, 0.1]
    filename_cat = 'truth.cat'
    _writeShears(n_angles,shears,filename_cat)
    logger.info('wrote file %s' % filename_cat)


    # this is the test config
    n_angles = 4
    shears = [0.05]
    filename_cat = 'truth.test.cat'
    _writeShears(n_angles,shears,filename_cat)
    logger.info('wrote file %s' % filename_cat)

def getGalaxyImages():

    config1 = copy.deepcopy(config)

    # get the total number of objects
    n_obj = len(galsim.config.GetNObjForMultiFits(config1,0,0))

    logger.info('found %d objects in the config file' % n_obj)

    # make sure that we are using correct obj_num and nimages
    if config1['args'].obj_num < 0: raise ValueError('start_id should be non-negative')
    else: obj_num = config1['args'].obj_num 
    if config1['args'].nimages <= 0: nimages = n_obj

    logger.info('starting with obj_num=%6d and processing %6d objects' % (obj_num,nimages))    
    
    # process imput  `
    galsim.config.ProcessInput(config1)

    logger.info('getting images')

    # get the galsim galaxy
    img_gal,img_psf,_,_ = galsim.config.BuildImages(config=config1,obj_num=obj_num,nimages=nimages,make_psf_image=True,logger=logger_config)
    logger.info('got %d images' % len(img_gal))

    if config['args'].verbosity > 2:

        import pylab
        # plot the all shear and all angles for 2 first galaxies
        for i in range(n_obj):
            pylab.subplot(1,2,1)
            pylab.imshow(img_gal[i].array,interpolation='nearest')
            pylab.subplot(1,2,2)
            pylab.imshow(img_psf[i].array,interpolation='nearest')
            filename_fig = 'debug/fig.buildimages.%03d.png' % i
            pylab.savefig(filename_fig)

    return img_gal

def runIm3shape():

    # open the RGC
    filepath_rgc = os.path.join(config['input']['real_catalog']['dir'],config['input']['real_catalog']['file_name'])
    rgc = pyfits.open(filepath_rgc)[1]

    # open the ring test catalog
    filename_cat = os.path.join(config['input']['catalog']['dir'],config['input']['catalog']['file_name'])
    ring_test_cat = numpy.loadtxt(filename_cat)

    # get im3shape
    dirpath_im3shape = os.path.join(os.environ['IM3SHAPE'],'python')
    sys.path.append(dirpath_im3shape)
    import im3shape

    # get images
    img_gals = getGalaxyImages()

    # get options
    n_pix = config['image']['size']
    i3o = im3shape.I3_options()
    i3o.read_ini_file(config['args'].filepath_ini)
    i3o.stamp_size = n_pix

    # get the file 
    filename_results = 'results.%s.%05d.cat' % (config['args'].filename_config,config['args'].obj_num)
    file_results = open(filename_results,'w')

    # create PSF from Moffat parameters
    # get i3 images - get first i3_galaxy to initialise the PSF - kind of strange, but hey..
    i3_galaxy = im3shape.I3_image(n_pix, n_pix)
    psf_beta = config['psf']['beta']
    psf_fwhm = config['psf']['fwhm']
    psf_e1 = config['psf']['ellip']['g1']
    psf_e2 = config['psf']['ellip']['g2']
    i3_psf = i3_galaxy.make_great10_psf(psf_beta, psf_fwhm, psf_e1, psf_e2, i3o)

    # loop over all created images
    for ig,img_gal in enumerate(img_gals):

        # get i3 images
        i3_galaxy = im3shape.I3_image(n_pix, n_pix)
        i3_galaxy.from_array(img_gal.array)  

        # get the galaxy id in the RGC
        id_ring = int(ring_test_cat[ig % config['settings']['n_repeat_gal'],columns['ring_test']['id_ring']])
        id_object = int(rgc.data['IDENT'][ig / config['settings']['n_repeat_gal']])
        unique_id = getUniqueID2(id_object,id_ring)

        i3_result, i3_best_fit = im3shape.i3_analyze(i3_galaxy, i3_psf, i3o, ID=unique_id)

        saveResult(file_results,i3_result)
        printResult(i3_result)

def saveResult(file_results,i3_result):

    fmt = '%d\t% e\t% 2.2f\t' + '% e\t'*13 + '%3d'*5 + '\n'
    line = fmt % (
                 i3_result.identifier,
                 i3_result.likelihood,
                 i3_result.time_taken,
                 i3_result.sersic_parameter_x0,
                 i3_result.sersic_parameter_y0,
                 i3_result.sersic_parameter_e1,
                 i3_result.sersic_parameter_e2,
                 i3_result.sersic_parameter_radius,
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
                 i3_result.levmar_residual_error_at_start,
                 i3_result.levmar_residual_error_at_end
                 )
    
    file_results.write(line)

def printResult(i3_result):

    fmt = '%d\t% e\t% 2.2f\t' + '% e\t'*6 
    line = fmt % (
                 i3_result.identifier,
                 i3_result.likelihood,
                 i3_result.time_taken,
                 i3_result.sersic_parameter_e1,
                 i3_result.sersic_parameter_e2,
                 i3_result.sersic_parameter_radius,
                 i3_result.sersic_bulge_flux,
                 i3_result.sersic_disc_flux,
                 i3_result.stats_signal_to_noise,
                 )
    
    logger.info(line)

if __name__ == "__main__":

    description = 'Noise and model bias driver. Requires $IM3SHAPE and $GALSIM to be set.'

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('command', type=str, help='command to run, available commands: {makering,run}')
    parser.add_argument('filepath_config', type=str, help='yaml config file, see reconvolution_validation.yaml for example.')
    parser.add_argument('--filepath_ini', type=str, default='nmb.ini', help='ini im3shape config file')
    parser.add_argument('--filepath_columns',  type=str, action='store', default='columns.yaml', help= 'columns file')
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-snr', '--signal_to_noise', type=float, action='store', default=1e20, help='signal to noise at which to run the test')
    parser.add_argument('--obj_num',  type=int, action='store', default= 0, help= 'first obj_num in config to process (starts from 1)') 
    parser.add_argument('--nimages',  type=int, action='store', default=-1, help= 'number of images to process, starting with obj_num')
    
    args = parser.parse_args()
    args.filename_config = os.path.basename(args.filepath_config)

    # Parse the integer verbosity level from the command line args into a logging_level string
    logging_levels = { 0: logging.CRITICAL, 
                       1: logging.WARNING,
                       2: logging.INFO,
                       3: logging.DEBUG }
    logging_level = logging_levels[args.verbosity]
    logging.basicConfig(format="%(message)s", level=logging_level, stream=sys.stdout)
    global logger , logger_config , columns , config
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
    config['gal']['signal_to_noise'] = args.signal_to_noise

    # load site config
    config['input']['real_catalog']['dir'] = os.path.join(os.environ['GALSIM'],'rgc')

    # load the columns file
    columns = yaml.load(open(args.filepath_columns,'r')) 

    # decide which command to run
    if args.command == 'makering':
        saveRingTestCatalog()
    elif args.command == 'run':
        saveRingTestCatalog()
        runIm3shape()
    else: raise ValueError('command %s not recognised')



