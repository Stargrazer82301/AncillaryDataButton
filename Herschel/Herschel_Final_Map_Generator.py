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
import gc
import datetime
import astropy.io.fits
import astropy.wcs
import astropy.io.votable
import astropy.modeling.models
import astropy.modeling.fitting
import aplpy
import ChrisFuncs
plt.ioff()





# Define function to finalise Herschel image of a given source in a given band
def Herschel_Generator(name, ra, dec, in_dir, band_dict, thumbnails):

    # Check if source has coverage in this band
    if not os.path.exists(in_dir+name+'_Herschel_'+band_dict['band']+'.fits'):
        print name+' has no coverage at '+band_dict['wavelength']
    else:

        # Read in image map
        image_in_fitsdata = astropy.io.fits.open(in_dir+name+'_Herschel_'+band_dict['band']+'.fits')
        image_in_cutout = image_in_fitsdata[0].data
        image_in_header = image_in_fitsdata[0].header
        image_in_wcs = astropy.wcs.WCS(image_in_header)
        image_out_cutout = image_in_cutout.copy()

        # Read in error map
        error_in_fitsdata = astropy.io.fits.open(in_dir+name+'_Herschel_'+band_dict['band']+'_Error.fits')
        error_in_cutout = error_in_fitsdata[0].data
        error_out_cutout = error_in_cutout.copy()

        # Set non-coverage pixels to be nan
        image_out_cutout[ np.where(image_out_cutout==0) ] = np.NaN
        error_out_cutout[ np.where(image_out_cutout==0) ] = np.NaN

        # Double-recheck that there is actual coverage at region of source
        i_target, j_target = np.round( image_in_wcs.wcs_world2pix( [[float(ra),float(dec)]], 0 )[0] ).astype(int).tolist()
        slice_target = image_out_cutout[ i_target-50:i_target+50, j_target-5:j_target+5 ]
        if np.where( np.isnan(slice_target)==False )[0].shape[0]==0:
            print name+' has no coverage at '+band_dict['wavelength']

        # If coverage present, continue with final processing
        else:

            # Convert SPIRE maps from MJy/sr to Jy/pix
            if band_dict['instrument'] == 'SPIRE':

                # Convert pixels from MJy/sr to Jy/sr
                image_out_cutout *= 1E6
                error_out_cutout *= 1E6

                # Calculate sr/pixel
                sr_per_sqarcsec = 2.3504E-11
                sr_per_pixels = sr_per_sqarcsec * band_dict['pix_size']**2

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
            out_header.set('TELESCOP', 'Herschel', 'Telescope that made this observation')
            out_header.set('INSTMNT', band_dict['instrument'], 'Instrument used for this observation')
            out_header.set('FILTER', band_dict['filter'], 'Filter used for this observation')
            out_header.set('WVLNGTH', band_dict['wavelength'], 'Wavelength of observation')
            out_header.set('DATABASE', 'Map produced using data acquired from the Herschel Science Archive (HSA)', 'Database from which the reduced data was acquired')
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
            image_cutout_hdulist.writeto(our_dir+name+'_Herschel_'+band_dict['band']+'.fits.gz', clobber=True)

            # Make thumbnail image of cutout
            if thumbnails:
                thumb_out = aplpy.FITSFigure(our_dir+name+'_Herschel_'+band_dict['band']+'.fits.gz')
                thumb_out.show_colorscale(cmap=band_dict['colourmap'], stretch='linear')#, pmin=7.5, pmax=99.5)
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.set_nan_color('black')
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=200, lw=2.0, edgecolor='#01DF3A', facecolor='#01DF3A' )
                thumb_out.save(our_dir+'Thumbnails/'+name+'_Herschel_'+band_dict['band']+'.png')
                thumb_out.close()

            # Close file, and clean memory
            image_in_fitsdata.close()
            error_in_fitsdata.close()
            gc.collect()





# Commence main task
if __name__ == "__main__":

    # State if thumbnail images are wanted
    thumbnails = False

    # State band information
    bands_dict = {'70':{'band':'70','instrument':'PACS','wavelength':'70um','filter':'PHOTBLUE','pix_size':2,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTBLUE','hdr_err_ext_name':'stDev'},
                  '100':{'band':'100','instrument':'PACS','wavelength':'100um','filter':'PHOTGREEN','pix_size':3,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTGREEN','hdr_err_ext_name':'stDev'},
                  '160':{'band':'160','instrument':'PACS','wavelength':'160um','filter':'PHOTRED','pix_size':4,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTRED','hdr_err_ext_name':'stDev'},
                  '250':{'band':'250','instrument':'SPIRE','wavelength':'250um','filter':'PSW','pix_size':6,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PSW','hdr_err_ext_name':'error'},
                  '350':{'band':'350','instrument':'SPIRE','wavelength':'350um','filter':'PMW','pix_size':8,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PMW','hdr_err_ext_name':'error'},
                  '500':{'band':'500','instrument':'SPIRE','wavelength':'500um','filter':'PLW','pix_size':12,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PLW','hdr_err_ext_name':'error'}}

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # Give paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Herschel/Mosaics/'
    our_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Herschel/Cutouts/'

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

        # In parallel, process each band
        print 'Processing '+str(name)
        pool = mp.Pool(processes=5)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( Herschel_Generator, args=(name, in_dir, band_dict, thumbnails,) )#
            Herschel_Generator(name, ra, dec, in_dir, band_dict, thumbnails)
        pool.close()
        pool.join()

    # Jubilate
    print 'All done!'
