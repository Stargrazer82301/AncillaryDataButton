# Identify location
import socket
location = socket.gethostname()
if location == 'saruman':
     dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Import smorgasbord
import os
import sys
import numpy as np
import scipy.ndimage
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





# Define function to provide appropriate Padding for FITS header entires
def Padding(entry):
    return ( ' ' * np.max([ 0, 19-len( str(entry) ) ]) ) + ' / '



# Define function to finalise Spitzer image of a given source in a given band
def Spitzer_Generator(name, ra, dec, in_dir, band_dict, thumbnails):
    band_short = band_dict['band_short']
    band_long = band_dict['band_long']
    wavelength = band_dict['wavelength']
    pix_size = band_dict['pix_size']
    colourmap = band_dict['colourmap']

    # Check if source has coverage in this band
    if not os.path.exists(in_dir+name+'_Spitzer_'+band_short+'.fits.gz'):
        print name+' has no coverage at '+wavelength
    else:

        # Read in image map
        image_in_fitsdata = astropy.io.fits.open(in_dir+name+'_Spitzer_'+band_short+'.fits.gz')
        image_in_cutout = image_in_fitsdata[0].data
        image_in_header = image_in_fitsdata[0].header
        image_in_wcs = astropy.wcs.WCS(image_in_header)
        image_out_cutout = image_in_cutout.copy()

        # Read in error map
        error_in_fitsdata = astropy.io.fits.open(in_dir+name+'_Spitzer_'+band_short+'_Error.fits.gz')
        error_in_cutout = error_in_fitsdata[0].data
        error_out_cutout = error_in_cutout.copy()

        # Set non-coverage pixels to be nan
        image_out_cutout[ np.where(image_out_cutout==0) ] = np.NaN
        error_out_cutout[ np.where(image_out_cutout==0) ] = np.NaN

        # Double-recheck that there is actual coverage at region of source
        i_target, j_target = np.round( image_in_wcs.wcs_world2pix( [[float(ra),float(dec)]], 0 )[0] ).astype(int).tolist()
        slice_target = image_out_cutout[ i_target-50:i_target+50, j_target-5:j_target+5 ]
        if np.where( np.isnan(slice_target)==False )[0].shape[0]==0:
            print name+' has no coverage at '+wavelength
        else:

            # Convert pixels from MJy/sr to Jy/sr
            image_out_cutout *= 1E6
            error_out_cutout *= 1E6

            # Calculate sr/pixel
            sr_per_sqarcsec = 2.3504E-11
            sr_per_pixels = sr_per_sqarcsec * pix_size**2

            # Convert pixels from Jy/sr to Jy/pix
            image_out_cutout *= sr_per_pixels
            error_out_cutout *= sr_per_pixels

            # Subtract minimum pixel value from entire map, to prevent issues with negative pixels later on
            image_out_cutout -= np.nanmin(image_out_cutout)

            # Only include coverage region that contains target source (ie, remove other, discontiguous regions)
            cont_binary = np.zeros([(image_out_cutout.shape)[0], (image_out_cutout.shape)[1]])
            cont_binary[ np.where( np.isnan(image_out_cutout)==False ) ] = 1
            cont_label = np.zeros([(image_out_cutout.shape)[0], (image_out_cutout.shape)[1]])
            cont_structure = np.array([[1,1,1], [1,1,1], [1,1,1]])
            scipy.ndimage.measurements.label(cont_binary, structure=cont_structure, output=cont_label)
            target_label = cont_label[i_target, j_target]
            image_out_cutout[ np.where( cont_label!=target_label ) ] = np.NaN
            error_out_cutout[ np.where( cont_label!=target_label ) ] = np.NaN

            # Create standard header
            out_header = astropy.io.fits.Header()
            date = datetime.datetime.now().isoformat()

            # Populate standard header entries
            out_header.set('TARGET', name , 'Target source of this map')
            out_header.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
            out_header.set('SIGUNIT', 'Jy/pix', 'Unit of the map')
            out_header.set('TELESCOP', 'Spitzer', 'Telescope that made this observation')
            out_header.set('INSTMNT', band_dict['instrument'], 'Instrument used for this observation')
            out_header.set('FILTER', band_long, 'Filter used for this observation')
            out_header.set('WVLNGTH', wavelength, 'Wavelength of observation')
            out_header.set('DATABASE', 'Map produced using data acquired from the Spitzer Heritage Archive (SHA)', 'Database from which the reduced data was acquired')
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
            image_cutout_hdulist.writeto(our_dir+name+'_Spitzer_'+band_short+'.fits.gz', clobber=True)

            # Make thumbnail image of cutout
            if thumbnails:
                thumb_out = aplpy.FITSFigure(our_dir+name+'_Spitzer_'+band_short+'.fits.gz')
                thumb_out.show_colorscale(cmap=colourmap, stretch='linear')#, pmin=7.5, pmax=99.5)
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.set_nan_color('black')
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=200, lw=2.0, edgecolor='#01DF3A', facecolor='#01DF3A' )
                thumb_out.save(our_dir+'Thumbnails/'+name+'_Spitzer_'+band_short+'.png')
                thumb_out.close()

            # Close file, and clean memory
            image_in_fitsdata.close()
            error_in_fitsdata.close()
            gc.collect()





# Commence main task
if __name__ == "__main__":

    # State if thumbnail images are wanted
    thumbnails = False

    # Decide what intrument to work for
    instrument = 'IRAC'

    # State band information
    if instrument=='IRAC':
        bands_dict = {'3.6um':{'instrument':'IRAC','band_short':'3.6','band_long':'IRAC 1','wavelength':'3.6um','pix_size':0.6,'colourmap':'pink'},
                      '4.5um':{'instrument':'IRAC','band_short':'4.5','band_long':'IRAC 2','wavelength':'4.5um','pix_size':0.6,'colourmap':'pink'},
                      '5.8um':{'instrument':'IRAC','band_short':'5.8','band_long':'IRAC 3','wavelength':'5.8um','pix_size':0.6,'colourmap':'pink'},
                      '8.0um':{'instrument':'IRAC','band_short':'8.0','band_long':'IRAC 4','wavelength':'8.0um','pix_size':0.6,'colourmap':'pink'}}
    elif instrument=='MIPS':
        bands_dict = {'24um':{'instrument':'MIPS','band_short':'24','band_long':'MIPS 1','wavelength':'24um','pix_size':2.45,'colourmap':'gist_heat'},
                      '70um':{'instrument':'MIPS','band_short':'70','band_long':'MIPS 2','wavelength':'70um','pix_size':4.0,'colourmap':'gist_heat'},
                      '160um':{'instrument':'MIPS','band_short':'160','band_long':'MIPS 3','wavelength':'160um','pix_size':8.0,'colourmap':'gist_heat'}}

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # Give paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Spitzer/Mosaics_'+instrument+'/'
    our_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Spitzer/Cutouts/'

    # Record time taken
    time_total = 0.0
    source_counter = 0.0
    source_total = name_list.shape[0]
    time_source_list = []

    # Loop over each source
    for i in range(0, name_list.shape[0]):#np.where(name_list=='NGC3031')[0]:
        source_counter += 1.0
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        name.replace('[','')
        name.replace(']','')

        # In parallel, process each band
        print 'Processing '+str(name)
        pool = mp.Pool(processes=5)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( Spitzer_Generator, args=(name, in_dir, band_dict, thumbnails,) )#
            Spitzer_Generator(name, ra, dec, in_dir, band_dict, thumbnails)
        pool.close()
        pool.join()

    # Jubilate
    print 'All done!'
