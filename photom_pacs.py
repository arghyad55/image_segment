import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
files=['hpacs_25HPPUNIMAPB_blue_1831_m0944_00_v1.0_1471689318484.fits','hpacs_25HPPUNIMAPR_1831_m0944_00_v1.0_1471689378571.fits']
cent_pix_x=[187.46,187.56]
cent_pix_y=[182.53,180.48]
beam_fwhm=[5.9,11.6]

for i in range(len(files)):
	image_file = get_pkg_data_filename(files[i])
	#print fits.info(image_file)
	image_data = fits.getdata(image_file,ext=1)
	header_data = fits.getheader(image_file,ext=1)
	data = image_data[650:1399,2000:2749]	

	from astropy.stats import SigmaClip
	from photutils.background import Background2D, MedianBackground
	sigma_clip = SigmaClip(sigma=3.)
	bkg_estimator = MedianBackground()
	bkg = Background2D(data, (30, 30), filter_size=(5,5), sigma_clip=sigma_clip,
         	          bkg_estimator=bkg_estimator)
	data -= bkg.background  # subtract the background
	threshold = 3. * bkg.background_rms  # above the background

	from astropy.convolution import Gaussian2DKernel, convolve
	from astropy.stats import gaussian_fwhm_to_sigma
	from photutils.segmentation import detect_sources, deblend_sources
	sigma = beam_fwhm[i] * gaussian_fwhm_to_sigma 
	kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
	convolved_data = convolve(data, kernel, normalize_kernel=True)
	npixels = 5
	segm = detect_sources(convolved_data, threshold, npixels=npixels)
	segm_deblend = deblend_sources(convolved_data, segm, npixels=npixels,
        	                       nlevels=32, contrast=0.001)

	from photutils import source_properties, EllipticalAperture
	from photutils.aperture import aperture_photometry
	labels=np.arange(63,68,1)
	cat = source_properties(data, segm,labels=labels)
	r = 3. # approximate isophotal extent
	apertures = []
	for obj in cat:
		position = (obj.xcentroid.value, obj.ycentroid.value)
		a = obj.semimajor_axis_sigma.value * r
		b = obj.semiminor_axis_sigma.value * r
		theta = obj.orientation.value
		apertures.append(EllipticalAperture(position, a, b, theta=theta))
	meas_pos = (cent_pix_x[i],cent_pix_y[i])
	meas_a=12.0639834607*r
	meas_b=9.56934286639*r
	meas_theta=1.07206495386
	meas_aperture= EllipticalAperture(meas_pos,meas_a,meas_b,theta=meas_theta)
	meas_tbl=aperture_photometry(data,meas_aperture)
	columns=['id','xcentroid','ycentroid','source_sum','area','semiminor_axis_sigma','semimajor_axis_sigma','orientation']
	tbl = cat.to_table(columns=columns)
	#print meas_tbl
	print tbl
	#tbl.write("pacs_aptr.html")

	plt.figure(figsize=(10, 8),dpi=240)
	plt.imshow(data,cmap='magma_r',origin='lower',norm=colors.PowerNorm(gamma=0.5))
	for aperture in apertures:
		aperture.plot(color='blue', lw=1.5, alpha=0.5)
	#meas_aperture.plot(color='#d62728')
	plt.colorbar()
	plt.clim(-0.01,5)
	bands=['pacs70u','pacs160u']
	#plt.savefig("/media/jdey/JYOTIRMOY/result_display/bgSubtracted_"+str(bands[files.index(items)])+".png")
	plt.show()
	

	'''header_data['CRPIX1']=1
	header_data['CRPIX2']=1
	header_data['CRVAL1']=277.620
	header_data['CRVAL2']=-10.911
	hdu = fits.PrimaryHDU(bkg.background,header=header_data)
	hdul_new = fits.HDUList([hdu])
	hdul_new.writeto('G20.989_pacs70_background.fits', overwrite = True)'''
