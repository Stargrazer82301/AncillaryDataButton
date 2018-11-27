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
bands = ['J','H','Ks']
wavelengths = ['1.25', '1.65', '2.16']

# Read in source catalogue
NESS_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
name_list = NESS_cat['name']
ra_list = NESS_cat['ra']
dec_list = NESS_cat['dec']

# Give paths
in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/2MASS/Mosaics/'
out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/2MASS/Cutouts/'

# Initiate timing
time_list = [ time.time() ]

# Loop over each source
for i in range(0, name_list.shape[0]):
    name = name_list[i].replace(' ','_')
    ra = ra_list[i]
    dec = dec_list[i]

    # Check if 2MASS data present
    nomap = 0
    for band in bands:
        if not os.path.exists(in_dir+str(name)+'_2MASS_'+str(band.replace('Ks','K'))+'.fits'):
            nomap += 1
    if nomap==3:
        continue
    time_start = time.time()

    # Loop over each band
    print 'Processing '+str(name)
    for b in range(0, len(bands)):
        band = bands[b]
        wavelength = wavelengths[b]

        # Check if data has already been processed; if so, skip
        if os.path.exists( os.path.join(out_dir, name+'_2MASS_'+band+'.fits' ) ):
            print 'Map for '+name+'_2MASS_'+band+'.fits already processed'
            continue

        # Check if input data exists exists; if not, skip
        if not os.path.exists( os.path.join( in_dir, name+'_2MASS_'+band.replace('Ks','K')+'.fits' ) ):
            print 'No data for '+name+'_2MASS_'+band+'.fits'
            continue

        # Read in map, skipping if file is corrupt
        try:
            in_image, in_header = astropy.io.fits.getdata( os.path.join( in_dir, name+'_2MASS_'+band.replace('Ks','K')+'.fits' ), header=True )
        except:
            print 'File '+name+'_2MASS_'+band.replace('Ks','K')+'.fits appears to be corrupted'
            continue
        in_wcs = astropy.wcs.WCS(in_header)
        out_image = in_image.copy()

        # Locate pixel coords of source
        in_wcs = astropy.wcs.WCS(in_header)
        location_pix = in_wcs.wcs_world2pix( np.array([[ np.float(ra), np.float(dec) ]]), 0 )[0]
        pix_i, pix_j = np.int(np.round(location_pix[1])), np.int(np.round(location_pix[0]))

        # Check map has coverage at source location; if not, continue
        if True in [ coord<=0 for coord in [ pix_i-5, pix_i+6, pix_j-5, pix_j+6 ] ]:
            continue
        try:
            in_image_slice = in_image[pix_i-5:pix_i+6, pix_j-5:pix_j+6]
        except:
            continue
        if np.where(np.isnan(in_image_slice)==True)[0].shape[0]==in_image_slice.size:
            continue

        # Create standard header
        cutout_header = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        cutout_header.set('TARGET', name , 'Target source of this map')
        cutout_header.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        cutout_header.set('SIGUNIT', 'Jy/pix', 'Unit of the map')
        cutout_header.set('TELESCOP', '2MASS', 'Telescope that made this observation')
        cutout_header.set('FILTER', band, 'Filter used for observation')
        cutout_header.set('WVLNGTH', wavelength+'um', 'Wavelength of observation')
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

        # Produce HDU
        image_cutout_hdu = astropy.io.fits.PrimaryHDU(data=out_image, header=image_cutout_header)

        # Create hdulist and save to file
        image_cutout_hdulist = astropy.io.fits.HDUList([image_cutout_hdu])
        image_cutout_hdulist.writeto( os.path.join(out_dir, name+'_2MASS_'+band+'.fits' ), clobber=True)

        # Make thumbnail image of cutout
        if thumbnails:
            thumb_out = aplpy.FITSFigure(out_dir+name+'_2MASS_'+band+'.fits')
            thumb_out.show_colorscale(cmap='gray')
            thumb_out.axis_labels.hide()
            thumb_out.tick_labels.hide()
            thumb_out.ticks.hide()
            thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=200, lw=2.0, edgecolor='#01DF3A', facecolor='#01DF3A' )
            thumb_out.save( os.path.join(out_dir, 'Thumbnails', name+'_2MASS_'+band+'.png') )
            thumb_out.close()

        # Close file, and clean memory
        gc.collect()

    # Report time until completion
    time_list.append( time.time() )
    time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
    print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'
