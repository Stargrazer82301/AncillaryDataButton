# Identify location
import socket
location = socket.gethostname()
if location == 'saruman':
    dropbox = '/home/user/spx7cjc/Desktop/Herdata/Dropbox/'

# Import smorgasbord
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import pdb
import gc
import time
import datetime
import shutil
import astropy.io.fits
import astropy.wcs
import astropy.io.votable
import astropy.modeling.models
import astropy.modeling.fitting
import aplpy
import ChrisFuncs
plt.ioff()



# Define function to provide appropriate Padding for FITS header entires
def Padding(entry):
    return ( ' ' * np.max([ 0, 19-len( str(entry) ) ]) ) + ' / '



# State if thumbnail images are wanted
thumbnails = False

# Provide band-specific information
bands = ['u','g','r','i','z']
wavelengths = ['355.10', '468.60', '616.60', '748.00', '893.20']

# Provide offsets between Vega and AB magnitudes (such that m_AB=m_SDSS+offset; hence m_SDSS=m_AB-offset)
SDSS_to_AB = [-0.04, 0.0, 0.0, 0.0, +0.02]

# Read in source catalogue
jingle_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
name_list = jingle_cat['name']
ra_list = jingle_cat['ra']
dec_list = jingle_cat['dec']

# Give paths
in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SDSS/Mosaics/'
out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SDSS/Cutouts/'

# Initiate timing
time_list = [ time.time() ]

# Loop over each source
for i in range(0, name_list.shape[0]):#np.where(name_list=='ESO359-002')[0]:
    name = name_list[i].replace(' ','_')
    ra = ra_list[i]
    dec = dec_list[i]

    # Check if SDSS data present
    nomap = 0
    for band in bands:
        if not os.path.exists(in_dir+str(name)+'_SDSS_'+str(band)+'.fits'):
            nomap += 1
    if nomap==5:
        continue
    time_start = time.time()

    # Loop over each band
    print 'Processing '+str(name)
    for b in range(0, len(bands)):
        band = bands[b]
        wavelength = wavelengths[b]
        """
        # Check if data has already been processed; if so, skip
        if os.path.exists( os.path.join(out_dir, name+'_SDSS_'+band+'.fits' ) ):
            print 'Map for '+name+'_SDSS_'+band+'.fits already processed'
            continue
        """
        # Read in map
        in_image, in_header = astropy.io.fits.getdata( in_dir+name+'_SDSS_'+band+'.fits', header=True )
        in_wcs = astropy.wcs.WCS( in_header )
        in_pix_width_arcsec = 3600.0 * ( np.min(np.abs(in_wcs.pixel_scale_matrix))**2.0 + np.max(np.abs(in_wcs.pixel_scale_matrix))**2.0 )**0.5
        out_image = in_image.copy()

        # Add minimum pixel value to image, to prevent NaNs from appearing later
        out_image = out_image + np.abs( 2.0 * np.nanmin( out_image ) )

        # Convert each pixel value to SDSS magnitudes
        out_image = 22.5 - ( 2.5 * np.log10(out_image) )

        # Convert each pixel value to AB magnitudes
        out_image = out_image + SDSS_to_AB[b]

        # Convert each pixel value to janskys
        out_image = ChrisFuncs.ABMagsToJy(out_image)

        # Scale pixel values to account pixel size increase from 0.396 to 0.45 arcseconds (as Montage conserves surface brightness)
        out_image *= in_pix_width_arcsec**2 / 0.396**2

        # Create standard header
        cutout_header = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        cutout_header.set('TARGET', name , 'Target source of this map')
        cutout_header.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        cutout_header.set('SIGUNIT', 'Jy/pix', 'Unit of the map')
        cutout_header.set('TELESCOP', 'SDSS', 'Telescope that made this observation')
        cutout_header.set('FILTER', band, 'Filter used for observation')
        cutout_header.set('WVLNGTH', wavelength+'nm', 'Wavelength of observation')
        cutout_header.set('MAPDATE', date, 'Date this map was produced by CJR Clark (Cardiff) from the reduced data')

        # Construct WCS system, and append to header
        cutout_wcs = astropy.wcs.WCS(naxis=2)
        cutout_wcs.wcs.crpix = [in_header['CRPIX1'], in_header['CRPIX2']]
        cutout_wcs.wcs.cdelt = [in_header['CDELT1'], in_header['CDELT2']]
        cutout_wcs.wcs.crval = [in_header['CRVAL1'], in_header['CRVAL2']]
        cutout_wcs.wcs.ctype = [in_header['CTYPE1'], in_header['CTYPE2']]
        cutout_wcs_header = cutout_wcs.to_header()
        for card in cutout_wcs_header.cards:
            cutout_header.append(card)

        # Produce header
        image_cutout_header = cutout_header.copy()

        # Save data to file
        image_cutout_hdulist = astropy.io.fits.writeto( out_dir+name+'_SDSS_'+band+'.fits', out_image, header=image_cutout_header, clobber=True )

        # Make thumbnail image of cutout
        if thumbnails:
            thumb_out = aplpy.FITSFigure(out_dir+name+'_SDSS_'+band+'.fits')
            thumb_out.show_colorscale(cmap='gist_gray')
            thumb_out.axis_labels.hide()
            thumb_out.tick_labels.hide()
            thumb_out.ticks.hide()
            thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=200, lw=2.0, edgecolor='#01DF3A', facecolor='#01DF3A' )
            thumb_out.save(out_dir+'Thumbnails/'+name+'_SDSS_'+band+'.png')
            thumb_out.close()

    # Report timings and clean memory
    time_list.append( time.time() )
    time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
    print 'Estimated completion time: '+time_est
    gc.collect()

# Jubilate
print 'All done!'




