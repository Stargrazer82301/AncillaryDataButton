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
bands = ['FUV','NUV']
wavelengths = ['152.8','227.1']

# Provide map zero-point magnitudes
zero_points = [18.82, 20.08]

# Read in source catalogue
jingle_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
name_list = jingle_cat['name']
ra_list = jingle_cat['ra']
dec_list = jingle_cat['dec']

# Give paths
in_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/GALEX/Mosaics/'
out_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/GALEX/Cutouts/'

# Initiate timing
time_list = [ time.time() ]

# Loop over each source
for i in range(0, name_list.shape[0]):#np.where(name_list=='NGC5584')[0]:
    name = name_list[i].replace(' ','_')
    ra = ra_list[i]
    dec = dec_list[i]

    # Check if GALEX data present
    nomap = 0
    for band in bands:
        if not os.path.exists(in_dir+str(name)+'_GALEX_'+str(band)+'.fits'):
            nomap += 1
    if nomap==2:
        continue

    # Loop over each band
    print 'Processing '+str(name)
    for b in range(0, len(bands)):
        band = bands[b]
        wavelength = wavelengths[b]
        zero_point = zero_points[b]

        # Check if GALEX data has already been processed; if so, skip
        if os.path.exists( os.path.join(out_dir, name+'_GALEX_'+band+'.fits' ) ):
            print 'Map for '+name+'_GALEX_'+band+'.fits already processed'
            continue

        # Check if input data exists exists; if not, skip
        if not os.path.exists( os.path.join( in_dir, name+'_GALEX_'+band+'.fits' ) ):
            print 'No data for '+name+'_GALEX_'+band+'.fits'
            continue

        # Read in map
        in_image, in_header = astropy.io.fits.getdata( in_dir+name+'_GALEX_'+band+'.fits', header=True )
        in_wcs = astropy.wcs.WCS(in_header)
        out_image = in_image.copy()

        # Convert each pixel value to AB magnitudes
        out_image = zero_point - ( 2.5 * np.log10(out_image) )

        # Convert each pixel value to janskys
        out_image = ChrisFuncs.ABMagsToJy(out_image)

        # Create standard header
        cutout_header = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        cutout_header.set('TARGET', name , 'Target source of this map')
        cutout_header.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        cutout_header.set('SIGUNIT', 'Jy/pix', 'Unit of the map')
        cutout_header.set('TELESCOP', 'GALEX', 'Telescope that made this observation')
        cutout_header.set('FILTER', band, 'Filter used for observation')
        cutout_header.set('WVLNGTH', wavelength+'nm', 'Wavelength of observation')
        cutout_header.set('MAPDATE', date, 'Date this map was produced by CJR Clark (Cardiff) from the reduced data')

        # Construct WCS system, and append to header
        cutout_wcs = astropy.wcs.WCS(naxis=2)
        cutout_wcs.wcs.crpix = [in_header['CRPIX1'], in_header['CRPIX2']]
        cutout_wcs.wcs.cdelt = [in_header['CD1_1'], in_header['CD2_2']]
        cutout_wcs.wcs.crval = [in_header['CRVAL1'], in_header['CRVAL2']]
        cutout_wcs.wcs.ctype = [in_header['CTYPE1'], in_header['CTYPE2']]
        cutout_wcs_header = cutout_wcs.to_header()
        for card in cutout_wcs_header.cards:
            cutout_header.append(card)

        # Produce header
        image_cutout_header = cutout_header.copy()

        # Save data to file
        image_cutout_hdulist = astropy.io.fits.writeto( out_dir+name+'_GALEX_'+band+'.fits', out_image, header=image_cutout_header, clobber=True )

        # Make thumbnail image of cutout
        if thumbnails:
            thumb_out = aplpy.FITSFigure(out_dir+name+'_GALEX_'+band+'.fits')
            thumb_out.show_colorscale(cmap='bone')
            thumb_out.axis_labels.hide()
            thumb_out.tick_labels.hide()
            thumb_out.ticks.hide()
            thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=200, lw=2.0, edgecolor='#01DF3A', facecolor='#01DF3A' )
            thumb_out.save(out_dir+'Thumbnails/'+name+'_GALEX_'+band+'.png')
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
