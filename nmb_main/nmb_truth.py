import os
import logging
import sys
import numpy
import sys
import math
import argparse
import yaml
import cPickle as pickle
import pyfits

dtype_table_truth   = { 'names'  : ['id_unique','id_cosmos','g1','g2','angle','id_angle','id_shear' , 'zphot'],
                        'formats': ['i8']*2 + ['f4']*3 + ['i4']*2 + ['f4']*1 }


def getUniqueID(id_object,id_angle,id_shear):
    return id_object*10000 + id_shear*100 + id_angle


def _writeShears(n_angles,shears,id_cosmos,redshift,file_cat):

    fmt = '%12d\t%8d\t% 2.4e\t% 2.4e\t% 2.10e\t%2d\t%2d\t%2.5f\n' 
    d_angle = 180./n_angles    
    for ig,g in enumerate(shears):
        shears_g12 = [ [ 0.0, +g ] , [ 0.0, -g ] , [ +g,   0.0 ] , [ -g,   0.0 ] , [ +g,  +g ] , [-g,  -g ] , [ +g,  -g ] , [ -g,  +g ]]
        for ig12,g12 in enumerate(shears_g12):
            for ia in range(n_angles):
                angle = ia*d_angle
                id_unique = getUniqueID(id_cosmos,ia,8*ig+ig12) 
                line = fmt % (id_unique, id_cosmos, g12[0], g12[1],   angle, ia, ig12, redshift); file_cat.write(line)

if __name__ == "__main__":

    description = 'Generate a truth table for nmb runs.'
    global logger,args

    # parse arguments
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('--filepath_ids',  type=str, action='store', help= 'filepath with ids of the galaxies to process')
    parser.add_argument('--filepath_out',  type=str, action='store', help= 'output file')
    parser.add_argument('--filepath_z', default='cosmos_zphot_mag25_id_z.tbl',   type=str, action='store', help= 'output file')
    parser.add_argument('--n_angles',  type=int, action='store', help= 'number of angles in the ring test')
    parser.add_argument('--shears',    type=float, nargs='+',  action='store', help= 'abs shear values, 8 will be created for each')
    args = parser.parse_args()

    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("nmb_truth") 

    logger.info('getting the truth table for %d angles and following shears' % args.n_angles)
    logger.info(str(args.shears))

    # load ids list
    ids_galsim = numpy.loadtxt(args.filepath_ids)
    data_cosmos = pyfits.open(args.filepath_z)[1].data
    ids_cosmos = data_cosmos['IDENT']
    ids_common = set(ids_cosmos).intersection(set(ids_galsim))
    logger.info('found %d common ids' % len(ids_common))

    # open the file
    file_cat = open(args.filepath_out,'w')
    # write header
    file_cat.write('# id_unique id_cosmos g1 g2 rotation_angle id_shear id_angle\n')
    # write the shears
    for i,idc in enumerate(ids_common):
        if (i % 100) == 0: logger.info('writing shears in galaxy %d' % i)
        redshift = data_cosmos[data_cosmos['IDENT']==idc]['ZPHOT']
        if redshift < 0: 
            logger.info('error redshift <0 %2.2f' % redshift)
        _writeShears(args.n_angles,args.shears,idc,redshift,file_cat) 
    file_cat.close()
    logger.info('finished writing file')
    truth_cat = numpy.loadtxt(args.filepath_out,dtype=dtype_table_truth) 
    n_obj = truth_cat.shape[0]
    logger.info('wrote file %s with %d objects' % (args.filepath_out,n_obj))

    # sort the truth cat
    truth_cat_sorted = truth_cat[numpy.argsort(truth_cat['id_unique'])]
    # write the pickled version
    filepath_pickle = args.filepath_out.replace('.cat','.pp')
    file_pickle = open(filepath_pickle,'w')
    pickle.dump(truth_cat_sorted,file_pickle,protocol=2)
    file_pickle.close()
    logger.info('wrote file %s with %d objects' % (filepath_pickle,n_obj))




