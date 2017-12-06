# Identify location
import socket
location = socket.gethostname()
if location == 'saruman':
    dropbox = '/home/user/spx7cjc/Desktop/Herdata/Dropbox/'

# Import smorgasbord
import os
import sys
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import pdb
import gc
import time
import datetime
import astropy.io.fits
import astropy.wcs
import astropy.io.votable
import astropy.modeling.models
import astropy.modeling.fitting
import aplpy
import ChrisFuncs
plt.ioff()



# State if thumbnail images are wanted
thumbnails = False

# Provide band-specific information
bands = ['W1','W2','W3','W4']
wavelengths = ['3.4', '4.6', '12', '22']

# Provide map zero-point magnitudes, and offsets between Vega and AB magnitudes (such that m_AB=m_WISE+offset; hence m_WISE=m_AB-offset)
zero_points = [20.5, 19.5, 18.0, 13.0]
WISE_to_AB = [+2.669, +3.339, +5.174, +6.620]

# Read in source catalogue
NESS_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
name_list = NESS_cat['name']
ra_list = NESS_cat['ra']
dec_list = NESS_cat['dec']

# Give paths
in_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/WISE/Mosaics/'
out_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/WISE/Cutouts/'

# Initiate timing
time_list = [ time.time() ]

# Loop over each source
for i in range(0, name_list.shape[0]):
    name = name_list[i].replace(' ','_')
    ra = ra_list[i]
    dec = dec_list[i]

    # Check if WISE data present
    nomap = 0
    for band in bands:
        if not os.path.exists(in_dir+str(name)+'_WISE_'+str(band)+'.fits'):
            nomap += 1
    if nomap==4:
        continue
    time_start = time.time()

    # Loop over each band
    print 'Processing '+str(name)
    for b in range(0, len(bands)):
        band = bands[b]
        wavelength = wavelengths[b]
        zero_point = zero_points[b]

        # Check if data has already been processed; if so, skip
        if os.path.exists( os.path.join(out_dir, name+'_WISE_'+wavelength+'.fits' ) ):
            print 'Map for '+name+'_WISE_'+band+'.fits already processed'
            continue

        # Read in image map
        image_in_fitsdata = astropy.io.fits.open(in_dir+name+'_WISE_'+band+'.fits')
        image_in_cutout = image_in_fitsdata[0].data
        image_in_header = image_in_fitsdata[0].header
        image_in_wcs = astropy.wcs.WCS(image_in_header)
        image_out_cutout = image_in_cutout.copy()

        # Read in error map
        error_in_fitsdata = astropy.io.fits.open(in_dir+name+'_WISE_'+band+'_Error.fits')
        error_in_cutout = error_in_fitsdata[0].data
        error_out_cutout = error_in_cutout.copy()

        # Add minimum pixel value to image map, to prevent NaNs from appearing later
        image_out_cutout = image_out_cutout + np.abs( 2.0 * np.nanmin( image_out_cutout ) )

        # Convert each pixel value to WISE Vega magnitudes
        image_out_cutout = zero_point - ( 2.5 * np.log10(image_out_cutout) )
        error_out_cutout = zero_point - ( 2.5 * np.log10(error_out_cutout) )

        # Convert each pixel value to AB magnitudes
        image_out_cutout = image_out_cutout + WISE_to_AB[b]
        error_out_cutout = error_out_cutout + WISE_to_AB[b]

        # Convert each pixel value to janskys
        image_out_cutout = ChrisFuncs.ABMagsToJy(image_out_cutout)
        error_out_cutout = ChrisFuncs.ABMagsToJy(error_out_cutout)

        # Create standard header
        out_header = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        out_header.set('TARGET', name , 'Target source of this map')
        out_header.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        out_header.set('SIGUNIT', 'Jy/pix', 'Unit of the map')
        out_header.set('TELESCOP', 'WISE', 'Telescope that made this observation')
        out_header.set('FILTER', band, 'Filter used for this observation')
        out_header.set('WVLNGTH', wavelength+'um', 'Wavelength of observation')
        out_header.set('DATABASE', 'Map produced using data acquired from the NASA/IPAC InfraRed Science Archive (IRSA)', 'Database from which the reduced data was acquired')
        out_header.set('MAPDATE', date, 'Date this map was produced by CJR Clark (Cardiff) from the reduced data')

        # Construct WCS for output header
        out_wcs = astropy.wcs.WCS(naxis=2)
        out_wcs.wcs.crpix = [image_in_header['CRPIX1'], image_in_header['CRPIX2']]
        out_wcs.wcs.cdelt = [image_in_header['CDELT1'], image_in_header['CDELT2']]
        out_wcs.wcs.crval = [image_in_header['CRVAL1'], image_in_header['CRVAL2']]
        out_wcs.wcs.ctype = [image_in_header['CTYPE1'], image_in_header['CTYPE2']]
        out_wcs_header = out_wcs.to_header()
        for card in out_wcs_header.cards:
            out_header.append(card)

        # Create image-specific header and HDU
        image_out_header = out_header.copy()
        image_out_header.set('EXTNAME', 'image')
        image_out_hdu = astropy.io.fits.PrimaryHDU(data=image_out_cutout, header=image_out_header)

        # Create error-specific header, with appropriate WCS
        error_out_header = out_header.copy()
        error_out_header.set('EXTNAME', 'error')
        error_out_hdu = astropy.io.fits.ImageHDU(data=error_out_cutout, header=error_out_header)

        # Create hdulist and save to file
        image_cutout_hdulist = astropy.io.fits.HDUList([image_out_hdu, error_out_hdu])
        image_cutout_hdulist.writeto(out_dir+name+'_WISE_'+wavelength+'.fits.gz', clobber=True)

        # Make thumbnail image of cutout
        if thumbnails:
            thumb_out = aplpy.FITSFigure(out_dir+name+'_WISE_'+wavelength+'.fits')
            thumb_out.show_colorscale(cmap='pink')
            thumb_out.axis_labels.hide()
            thumb_out.tick_labels.hide()
            thumb_out.ticks.hide()
            thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=200, lw=2.0, edgecolor='#01DF3A', facecolor='#01DF3A' )
            thumb_out.save(out_dir+'Thumbnails/'+name+'_WISE_'+wavelength+'.png')
            thumb_out.close()

   # Report timings and clean memory
    time_list.append( time.time() )
    time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
    print 'Estimated completion time: '+time_est
    gc.collect()

# Jubilate
print 'All done!'





"""
# Define median filter function
def clip_filter(image_in, kernel):
    image_ds = downsample(image_in, 5)
    mirror = mirror_grid(image_ds)
    kernel = np.int( np.round( ( np.float(kernel) / 10.0 ) ) )
    kernel = np.int( np.round(2.0*kernel) ) - np.int( np.abs( np.mod( np.round(2.0*kernel), 2.0 ) - 1.0 ) )
    image_out = np.zeros([image_ds.shape[0], image_ds.shape[1]])
    for i in range(0, image_ds.shape[0]):
        print 'Row '+str(i+1)+' of '+str(image_ds.shape[0])
        for j in range(0, image_ds.shape[1]):
            image_out[i,j] = slicing_clip(image_ds, mirror, kernel, i, j)[0]
    image_out = congrid.congrid(image_out, [image_in.shape[0],image_in.shape[1]], minusone=True)
    return image_out
"""

"""
# Define median filter function
def clip_filter(image_in, kernel):
    image_ds = downsample(image_in, 5)
    mirror = mirror_grid(image_ds)
    kernel = np.int( np.round( ( np.float(kernel) / 10.0 ) ) )
    kernel = np.int( np.round(2.0*kernel) ) - np.int( np.abs( np.mod( np.round(2.0*kernel), 2.0 ) - 1.0 ) )
    list_out = [ [ np.nan ]*image_ds.shape[1] ]*image_ds.shape[0]
    pool = mp.Pool(processes=50)
    for i in range(0, image_ds.shape[0]):
        for j in range(0, image_ds.shape[1]):
            list_out[i][j] = pool.apply_async( slicing_clip, args=(image_ds, mirror, kernel, i, j,) )
    pool.close()
    pool.join()
    pdb.set_trace()
    return list_out
"""
