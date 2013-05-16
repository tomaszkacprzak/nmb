import os
from numpy import *
from StringIO import StringIO 
import copy
import matplotlib.pyplot as plt
import pyfits
import subprocess

# dir_im3 = os.path.join(os.path.split(__file__)[0], "..")
dir_im3 = '/home/kacprzak/code/ucl_des_shear/im3shape/'

dir_out = '.'
file_img = "temp.img.fits"
cat_gal  = "temp.gal.cat"
cat_psf  = "temp.psf.cat"
cat_im3  = "temp.im3.cat"
conf_ini = "temp.ini"
file_img_psf = "temp.psf.fits"



class Galaxy:
	def __init__(self):
		self.model_gal = "sersics"
		self.params_gal_true = array([19.5, 19.5, 0.3, 0.2, 2., 1., 0.8, 0.2, 4., 1., 0., 0.])
		self.params_gal_measured = ones([12])*66.
		self.params_all_measured = ones([20])*66.
		self.model_psf = "moffat"
		self.params_psf_true = array([3., 2.85 ,0., 0.])
		self.snr = 20
		self.image_true = zeros([39,39])
		self.image_noisy = zeros([39,39])
		self.image_bestfit = zeros([39,39])
		self.image_psf_hires = "none"
		self.params_gal_labels=("x0","y0","e1","e2","r","r_ratio","a1","a2","index_bulge","index_disc","delta_e","delta_theta")
		self.load_ini_file()
		
	def get_image_true(self,file_ini="default",file_fits=file_img,show_cmd=False):
		snr = self.snr
		self.snr = -1
		image = self.get_image(file_ini,file_fits,show_cmd)
		self.snr = snr
		self.image_true = image
		return self.image_true
		
	def get_image_noisy(self,file_ini="default",file_fits=file_img,show_cmd=False):
		image = self.get_image(file_ini,file_fits,show_cmd=False)
		self.image_noisy = image
		return self.image_noisy
	
	def get_image(self,file_ini="default",file_fits=file_img,show_cmd=False):
		
		if file_ini=="default":
			
			file_ini = conf_ini
		else:
			self.load_ini_file()
			
		self.write_temp_files()
			
		if self.image_psf_hires == "none":
				
			params_gal_str =  ''.join('%+2.4e ' % n for n in self.params_gal_true)
			params_psf_str =  ''.join('%+2.4e ' % n for n in self.params_psf_true)
			
			#image = loadtxt("%s/%s" % (dir_out,file_fits))
		else:
						
			params_gal_str =  ''.join('%+2.4e ' % n for n in self.params_gal_true)
			params_psf_str = "file %s" % file_img_psf
			 
		cmd = "%s/bin/im3generate %s %s %s %s %s %s %s" % (dir_im3,file_ini,file_fits,str(self.snr),self.model_psf,params_psf_str,self.model_gal,params_gal_str)
		if show_cmd==True: print cmd
		
		# os.system(cmd)

		process = subprocess.Popen(cmd.split(' '))
		process.wait()

		image = pyfits.getdata("%s/%s" % (dir_out,file_fits))
			
				
		
		return image

	def get_measurement(self,file_ini="default",cg=cat_gal,cp=cat_psf,cm=cat_im3,fi=file_img,fip=file_img_psf,show_cmd=False):
		
		if file_ini=="default":
			file_ini = conf_ini
		else:
			self.load_ini_file()		
			
		self.write_temp_files()
						
		if self.image_psf_hires == "none":
						
			cmd = "%s/bin/im3shape %s %s %s %s %s" % (dir_im3,file_ini,fi,cg,cp,cm)
			if show_cmd==True: print cmd
			# os.system(cmd)

			process = subprocess.Popen(cmd.split(' '))
			process.wait()
				
		else:
			cmd = "%s/bin/im3shape %s %s %s %s %s" % (dir_im3,file_ini,fi,cg,fip,cm)
			if show_cmd==True: print cmd
			# os.system(cmd)
	
			process = subprocess.Popen(cmd.split(' '))
			process.wait()
		
		self.params_all_measured = loadtxt(cm)
		self.params_gal_measured = copy.deepcopy(self.params_gal_true)
		self.params_gal_measured[0:5] = self.params_all_measured[2:7]
		self.params_gal_measured[6:8] = self.params_all_measured[7:9]
		
		if os.path.isfile('model.1.fits'):
			self.image_bestfit =  pyfits.getdata("model.1.fits")
		else:
			print 'bestfit image not available'

	def write_temp_files(self,cg=cat_gal,cp=cat_psf,cm=cat_im3,fi=file_img,fip=file_img_psf,ini=conf_ini):
		
		# write the i3_catalog
		f = open(cg,'w')
		f.write("1 %2.4f %2.4f" % (self.params_gal_true[0],self.params_gal_true[1]))
		f.close()

		# write the i3_psf
		f = open(cp,'w')
		f.write("%2.4f %2.4f %2.4f %2.4f " % (self.params_psf_true[0],self.params_psf_true[1],self.params_psf_true[2],self.params_psf_true[3]))
		f.close()
		
		# write the ini file		
		f = open(ini,'w')
		f.write(self.ini)
		f.close()
		
		# write the noisy image
		hdu = pyfits.PrimaryHDU(self.image_noisy)
		os.system('rm ' + fi)
		hdu.writeto(fi)
		
		if self.image_psf_hires != "none":
		
			hdu = pyfits.PrimaryHDU(self.image_psf_hires)
			os.system('rm ' + fip)
			hdu.writeto(fip)
		
		# savetxt(fi,self.image_noisy)
					
	def load_ini_file(self,file_ini="default"):
	
		if file_ini=="default":	
			self.ini = "model_name = sersics\nstamp_size = 39\nupsampling = 3\npadding = 2\nverbosity = -1\nminimizer_verbosity = -1\nminimizer_max_iterations=500\nlevmar_eps1=1e-40\nlevmar_eps2=1e-10\nlevmar_eps3=1e-40\nlevmar_eps4=1e-10\nlevmar_tau=1e-4\nlevmar_LM_INIT_MU=1e-4\nn_central_pixels_to_upsample= 6\nn_central_pixel_upsampling= 15\nsave_images=yes\nbackground_subtract=no\nrescale_stamp=no\npsf_truncation_pixels=10000\nminimizer_loops=3\n"
		else: 
			f = open(file_ini,'r')
			self.ini = f.read()
			f.close()
		
	def get_fwhm(self):
		
		fwxm=0.5
		n_pix=39
		n_sub=9
		image_nx = n_pix*n_sub
		xy0 = n_pix/2 + 0.5;
					
		gal.ini = "model_name = sersics\nstamp_size = 39\nupsampling = 9\npadding = 0\nverbosity = -1\nminimizer_verbosity = -1\nminimizer_max_iterations=500\nlevmar_eps1=1e-40\nlevmar_eps2=1e-10\nlevmar_eps3=1e-40\nlevmar_eps4=1e-10\nlevmar_tau=1e-4\nlevmar_LM_INIT_MU=1e-4\nn_central_pixels_to_upsample= 6\nn_central_pixel_upsampling= 15\nsave_images=yes\nbackground_subtract=no\nrescale_stamp=no\npsf_truncation_pixels=10000\nminimizer_loops=3\nimage_generator_resolution_high=yes\n"
		gal.params_gal_true[0]=xy0
		gal.params_gal_true[1]=xy0
		gal.params_gal_true[2]=0.
		gal.params_gal_true[3]=0.
		
		gal.get_image_true()
		
		gal.image_true = gal.image_true - amin(gal.image_true.flatten()) 
				
		f1 = 0.
		f2 = 0.
		x1 = 0
		x2 = 0
		
				
		profile = gal.image_true[int(image_nx//2),:]
				
		max_ind = int(image_nx//2)
		max_val = profile[max_ind]
		f3 = max_val*fwxm
		
		diff = abs(profile-f3)
		
		x1 = argmin(diff)
		f1 = profile[x1]
		
		if( f1 < f3 ): 	x2 = x1+1
		else: 		x2 = x1-1
		f2 = profile[x2];
	
		a = (f1-f2)/(x1 - x2);
		b = f1 - a*x1; 
		x3 = (f3 - b)/a;
		
		
					
		fwhm = 2.*abs(max_ind-x3);
		#fwhm = 2.*abs(xy0*n_sub-x3);
		Rgp = fwhm/n_sub
				
		print (f1,f2,f3,x1,x2,x3,n_sub,xy0,xy0*n_sub,max_ind,argmax(profile),fwhm,Rgp)
		
		#if isinf(Rgp) or isnan(Rgp): 
		plt.subplot(211)
		plt.imshow(gal.image_true,interpolation="nearest")
		plt.subplot(212)
		plt.plot(arange(image_nx)/float(n_sub),profile)
		plt.plot(float(x1)/float(n_sub),f1,'rx')
		plt.plot(float(x2)/float(n_sub),f2,'gx')
		plt.plot(float(x3)/float(n_sub),f3,'bx')
		plt.plot(float(max_ind)/float(n_sub),max_val,'bx')
		plt.show()
		
		#Rgp=666.
		
		return Rgp
		
	def write_image(self,file_fits, image):
		hdu = pyfits.PrimaryHDU(image)
		hdu.writeto(file_fits,clobber=True);

	def write_image_true(self,file_fits):
		self.write_image(file_fits,self.image_true)

	def write_image_noisy(self,file_fits):
		self.write_image(file_fits,self.image_noisy)
		
	def get_string_from_vector(self,vector,fmt='%2.4e'):
		string = ''
		for i in vector:
			string += fmt % i
			string += '\t'
		string+='\n'
		return string


#class Catalog:
	#def __init__(self,filename_catalog,filename_ini="default"):
		
		#self.filename_cat = filename_catalog
		#self.filename_ini = filename_ini
		#self.list_gals = self.load_catalog(filename_catalog,filename_ini)
				
def load_catalog(file_cat,file_ini="default"):
	
	index_gal_params = arange(6,18)
	index_psf_params = arange(2,6)
	index_snr = 1
	
	print file_cat
	cat_sky_t = loadtxt(file_cat,ndmin=2)



	if cat_sky_t.shape[1] == 1:
		cat_sky = cat_sky_t.transpose()
	else: 	
		cat_sky = cat_sky_t
			
	n_gals = cat_sky.shape[0]
			
	print "found %d objects in catalogue %s" % (n_gals,file_cat)
		
	list_gals = [Galaxy() for ng in range(n_gals)];
	
	#for gal in list_gals:
	for ng in range(n_gals):
	
		gal = list_gals[ng]
		gal.params_gal_true = cat_sky[ng,:].take(index_gal_params)
		gal.params_psf_true = cat_sky[ng,:].take(index_psf_params)
		gal.snr = cat_sky[ng,index_snr]
		gal.load_ini_file(file_ini)
		
	#for ng in range(n_gals): print list_gals[ng].params_gal_true
		
	return list_gals
	
		
if __name__=="__main__":
	
	
	gal = Galaxy()
	gal.snr = -1
	gal.params_gal_true = array([19.5, 19.5, 0., 0., 2.7, 1., 0.5, 0.5, 4., 1., 0., 0.])
	#gal.params_psf_true[0] = 3.5
	
	gal.get_image_true(show_cmd=True)
	gal.get_image_noisy(show_cmd=True)
	
	plt.imshow(gal.image_noisy,interpolation="nearest")
	plt.show()
	plt.imshow(gal.image_true,interpolation="nearest")
	plt.show()
	
	#gal.get_measurement(show_cmd=True)
	print gal.get_fwhm()/2.85