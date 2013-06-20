import logging

def loadTable(table_name,filepath,dtype=None,hdu=1,logger=logging):

    try:
        table = eval(table_name)
    except:

        logger.debug('loading %s' % filepath)
        if filepath.split('.')[-1] == 'pp':
                import cPickle as pickle
                file_pickle = open(filepath)
                table = pickle.load(file_pickle)
                file_pickle.close()

        elif filepath.split('.')[-1] == 'fits' or filepath.split('.')[-2] == 'fits':
                import pyfits
                fits = pyfits.open(filepath)
                table = fits[hdu].data

                # if using FITS_record, get FITS rec
                if isinstance(table,pyfits.FITS_record):
                    import pdb; pdb.set_trace()
                    table = table.array
        else:
                import numpy
                table = numpy.loadtxt(filepath,dtype=dtype)
        
    else:
        logger.debug('using preloaded array %s' % table_name)
    
    logger.debug('loaded %s correctly, got %d rows' % (filepath,len(table)))

    globals()[table_name] = table
    return table

def saveTable(filepath,table,logger=logging):

    import numpy
    formats = { numpy.dtype('int64') : '% 12d' ,
                numpy.dtype('int32') : '% 12d' ,
                numpy.dtype('float32') : '% .10e' ,
                numpy.dtype('float64') : '% .10e' }

    if filepath.split('.')[-1] == 'pp':
        import cPickle as pickle
        file_pickle = open(filepath,'w')
        pickle.dump(table,file_pickle,protocol=2)
        file_pickle.close()
    elif filepath.split('.')[-1] == 'fits' or filepath.split('.')[-2] == 'fits':
        import pyfits
        if type(table) is pyfits.core.HDUList:
            table.writeto(filepath,clobber=True)
        else:
            fits = getFITSTable(table)
            fits.writeto(filepath,clobber=True)                      
    else:
        header = '# ' + ' '.join(table.dtype.names)
        fmt = [formats[table.dtype.fields[f][0]] for f in table.dtype.names]
        float(numpy.__version__[0:3])
        if float(numpy.__version__[0:3]) >= 1.7:
            numpy.savetxt(filepath,table,header=header,fmt=fmt,delimiter='\t')
        else:
            numpy.savetxt(filepath,table,fmt=fmt,delimiter='\t')
            with file(filepath, 'r') as original: data = original.read()
            with file(filepath, 'w') as modified: 
                modified.write(header + '\n' + data)
                modified.close()
            

    logger.info('table saved %s correctly, got %d rows' % (filepath,len(table.shape)))

def getFITSTable(numpy_array):

    import numpy
    import pyfits
    formats = { numpy.dtype('int64') : 'K' , numpy.dtype('int16') : 'I' , numpy.dtype('int32') : 'J' , numpy.dtype('float32') : 'E' , numpy.dtype('float64') : 'D' ,
                    numpy.dtype('>i8') : 'K', numpy.dtype('>i4') : 'I', numpy.dtype('>f4') : 'E' , numpy.dtype('>f8') : 'D'}

    cols = []

    for i,col_name in enumerate(numpy_array.dtype.names):
        col_type = numpy_array.dtype[i]
        col_fmt = formats[col_type]
        col = pyfits.Column(name=col_name,format=col_fmt,array=numpy_array[col_name])
        cols.append(col)

    tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
    hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([hdu, tbhdu])

    return hdulist

def getColNames(table):
    
    import numpy
    import pyfits

    if isinstance(table,pyfits.FITS_record):
        colnames = table.array.dtype.names
    elif isinstance(table,pyfits.FITS_rec):
        colnames = table.dtype.names
    elif isinstance(table,numpy.ndarray):
        colnames = table.dtype.names

    return colnames