dtype_table_truth   = { 'names'  : ['id_unique','id_cosmos','g1','g2','angle','id_angle','id_shear' , 'zphot'],
                        'formats': ['i8']*2   + ['f4']*3 + ['i4']*2 + ['f4']*1 }

dtype_table_results = { 'names'   : ['identifier','likelihood','time_taken','x0','y0','e1','e2','radius','fwhm','bulge_flux','disc_flux','flux_ratio','signal_to_noise','min_residuals','max_residuals','model_min','model_max','number_of_likelihood_evals','number_of_iterations','reason_of_termination'],
                        'formats' : ['i8'] + ['f4']*16 + ['i4']*3 }           

dtype_table_stats =  { 'names'   : ['index', 'cosmos_id' , 'zphot' ,'m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std' , 'hlr' , 'rgp' , 'snr'],
                        'formats' : ['i8']*2 + ['f4']*12 } 


dtype_table_results2 = { 'names'   : ['id_global' , 'id_object' , 'id_unique', 'id_cosmos', 'likelihood','time_taken','x0','y0','e1','e2','radius','fwhm','bulge_flux','disc_flux','flux_ratio','signal_to_noise','min_residuals','max_residuals','model_min','model_max','number_of_likelihood_evals','number_of_iterations','reason_of_termination'],
                        'formats' : ['i8']*4 + ['f4']*16 + ['i4']*3 }           

dtype_table_binstats = { 'names'   : ['bin_id', 'bin_value' , 'bin_ngals' ,'m1','m2','m1_std','m2_std','c1','c2','c1_std','c2_std'],
                        'formats' : ['i8']*1 + ['f4'] + ['i8'] + ['f4']*8 } 



NO_RESULT_FLAG = 666
