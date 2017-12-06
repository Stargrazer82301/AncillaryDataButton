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
def Planck_Generator(name, in_dir, band_dict, thumbnails):
    wavelength = band_dict['wavelength']

    # Check if source has coverage in this band
    if not os.path.exists(in_dir+name+'_Planck_'+wavelength+'.fits'):
        print name+' has no coverage at '+wavelength
    else:

        # Read in map
        in_fitsdata = astropy.io.fits.open(in_dir+name+'_Planck_'+wavelength+'.fits')
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_wcs = astropy.wcs.WCS(in_header)
        pix_width_arcsec = 3600.0 * np.abs( np.max( in_wcs.pixel_scale_matrix ) )
        out_image = in_image.copy()

        # Set non-coverage pixels to be nan
        out_image[ np.where(out_image==0) ] = np.NaN

        # Calculate sr/pixel
        sr_per_sqarcsec = 2.3504E-11
        sr_per_pixels = sr_per_sqarcsec * pix_width_arcsec**2

        # Convert pixel units from K(cmb) to MJy/sr
        if band_dict['units']=='K(cmb)':
            out_image *= band_dict['conversion']

        # Convert pixel units from MJy/sr to Jy/pix
        out_image *= 1E6
        out_image *= sr_per_pixels

        # Create standard header
        cutout_header = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        cutout_header.set('TARGET', name, 'Target source of this map')
        cutout_header.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        cutout_header.set('SIGUNIT', 'Jy/pix', 'Unit of the map')
        cutout_header.set('TELESCOP', 'Planck', 'Telescope that made this observation')
        cutout_header.set('INSTMNT', band_dict['instrument'], 'Instrument used for this observation')
        cutout_header.set('FREQ', band_dict['freq'], 'Filter used for this observation')
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
        image_cutout_hdulist.writeto(our_dir+name+'_Planck_'+wavelength+'.fits', clobber=True)

        # Make thumbnail image of cutout
        if thumbnails:
            try:
                thumb_out = aplpy.FITSFigure(our_dir+name+'_Planck_'+wavelength+'.fits')
                thumb_out.show_colorscale(cmap='gist_heat', stretch='arcsinh')
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=500, lw=2.5, edgecolor='#01DF3A')
                thumb_out.save(our_dir+'Thumbnails/'+name+'_Planck_'+wavelength+'.png')
                thumb_out.close()
            except:
                print 'Failed making thumbnail for '+name
                pdb.set_trace()

        # Close file, and clean memory
        in_fitsdata.close()
        gc.collect()





# Commence main task
if __name__ == "__main__":

    # State if thumbnail images are wanted
    thumbnails = False

    # State band information
    bands_dict = {'Planck030':{'band_name':'Planck 030','freq':'30GHz','wavelength':'10600','beam_area':1188.945,'units':'K(cmb)','instrument':'LFI','conversion':0.979328*24.845597},
                  'Planck044':{'band_name':'Planck 044','freq':'44GHz','wavelength':'6810','beam_area':832.168,'units':'K(cmb)','instrument':'LFI','conversion':0.95121302*59.666236},
                  'Planck070':{'band_name':'Planck 070','freq':'70GHz','wavelength':'4260','beam_area':200.519,'units':'K(cmb)','instrument':'LFI','conversion':0.88140690*151.73238},
                  'Planck100':{'band_name':'Planck 100','freq':'100GHz','wavelength':'3000','beam_area':105.777,'units':'K(cmb)','instrument':'HFI','conversion':0.76581996*306.81118},
                  'Planck143':{'band_name':'Planck 143','freq':'143GHz','wavelength':'2100','beam_area':59.952,'units':'K(cmb)','instrument':'HFI','conversion':0.59714682*627.39818},
                  'Planck217':{'band_name':'Planck 217','freq':'217GHz','wavelength':'1380','beam_area':28.426,'units':'K(cmb)','instrument':'HFI','conversion':0.31573332*1444.7432},
                  'Planck353':{'band_name':'Planck 353','freq':'353GHz','wavelength':'850','beam_area':26.653,'units':'K(cmb)','instrument':'HFI','conversion':0.071041398*3823.1434},
                  'Planck545':{'band_name':'Planck 545','freq':'545GHz','wavelength':'550','beam_area':26.305,'units':'MJy/sr','instrument':'HFI'},
                  'Planck857':{'band_name':'Planck 857','freq':'857GHz','wavelength':'350','beam_area':23.985,'units':'MJy/sr','instrument':'HFI'}}

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
    ness_cat = ness_cat[::-1]
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # Give paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Planck/Raw/'
    our_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Planck/Cutouts/'

    # Record time taken
    time_list = [time.time()]

    # Loop over each source
    for i in range(0, name_list.shape[0]):#np.where(name_list=='NGC3031')[0]:
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        name.replace('[','')
        name.replace(']','')

        # Check if Spitzer data present
        nomap = 0
        for band in bands_dict.keys():
            if not os.path.exists(in_dir+str(name)+'_Planck_'+str(band)+'.fits'):
                nomap += 1
        if nomap==5:
            continue

        # In parallel, process each band
        print 'Processing '+str(name)
        pool = mp.Pool(processes=9)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( Planck_Generator, args=(name, in_dir, band_dict, thumbnails,) )
            Planck_Generator(name, in_dir, band_dict, thumbnails)
        pool.close()
        pool.join()

        # Report predicted time remaining
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        print 'Estimated completion time: '+time_est

    # Jubilate
    print 'All done!'
