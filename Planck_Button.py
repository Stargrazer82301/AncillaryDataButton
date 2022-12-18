# Import smorgasbord
import os
import multiprocessing as mp
import numpy as np
import astropy.io.votable
import astropy.coordinates
from astroquery.skyview import SkyView
import aplpy
import urllib
import gc
import re
import shutil
import copy
import pdb
import time
import datetime
import ChrisFuncs





# Define main function
def Run(ra, dec, width, name=None, out_dir=None, temp_dir=None, replace=False, flux=True, thumbnails=False):
    """
    Function to generate standardised cutouts of Planck observations.

    Arguments
        ra: {float, sequence of float}
                A sequence of right ascension values, in decimal degrees, of the targets to be processed. Alternatively,
                if you're only interested in one target, a single RA value can be given here.
        dec: {float, sequence of float}
                A sequence of declination values, in decimal degrees, of the targets to be processed. Alternatively, if
                you're only interested in one target, a single Dec value can be given here.
        width: {float, sequence of float}
                A sequence giving the desired width of the cutout square for each target, in decimal degrees.
                Alternatively, if you're only interested in one target, a single width value can be given here.

    Keyword arguments
        name: {str, sequence of str}, optional
                A sequence giving the name of each target; if you're only interested in one target, a
                single name can be given here. If not provided, a name is constructed automatrically from the target
                coordinates, according to the IAU catalogue convention.
        out_dir: str, optional
                A string giving the path to the directory where the output FITS files will be placed. If not provided,
                files will simply be written to the current working directory.
        temp_dir: str, optional
                A string giving the path to be used as a temporary working directory by Planck_Button. If not provided,
                a temporary directory will be created inside the output directory.
        replace: bool, optional
                If False, Planck_Button will search the output directory for any pre-existing output FITS files from
                previous runs of the function, and will not bother repeat creating these maps (making it easy to resume
                processing a large number of targets from an interruption. If True, Planck_Button will produce maps for
                all input targets, regardless of whether maps for these targets already exist in the output directory.
        flux: bool, optional
                If True, output maps will be in flux density units of Jy/pix. If false, output maps will be in surface
                brightness units of MJy/sr.
        thumbnails: bool, optional
                If True, JPG thumbnail images of the generated maps will also be proced and placed in out_dir.
    """


    # Make sure input values are in list format, and sort out variable names for rest of function
    if not hasattr(ra, '__iter__'):
        ra = [ra]
    ra_list = np.array(ra)
    del(ra)
    if not hasattr(dec, '__iter__'):
        dec = [dec]
    dec_list = np.array(dec)
    del(dec)

    # Check that ra and declists all have same lengths
    if np.std([float(len(ra_list)), float(len(dec_list))]) > 0:
        raise Exception('Input sequences of ra and dec all need to be the same length')

    # If single width provided, but multiple coordinates, create width array of same value repeated required number of times
    if not hasattr(width, '__iter__'):
        if len(ra_list) > 1:
            width_list = [width] * len(ra_list)

        # Else, if only one RA and one width given, stick width value into list, too
        elif len(ra_list) == 1:
            width_list = [width]
    width_list = np.array(width_list)
    del(width)

    # If no names provided, use coordinates to generate standardised names as per IAU catalogue convention
    if not hasattr(name, '__iter__'):
        if (name == None):
            name = []
            for i in range(len(ra_list)):
                coord = astropy.coordinates.SkyCoord(str(ra_list[i])+'d '+str(dec_list[i])+'d')
                name_coord = re.sub('[hmsdms. ]', ' ', coord.to_string('hmsdms'))
                name_coord = name_coord.split(' ')
                name_coord[3] = name_coord[3][:min(2,len(name_coord[3]))]
                name_coord[8] = name_coord[8][:min(2,len(name_coord[8]))]
                name_coord = 'J'+''.join(name_coord)
                name.append(re.sub('[hmsdms. ]', ' ', coord.to_string('hmsdms')))

        # If only one name provided, stick it into an array
        name_list = np.array([name])

    # If a sequence of names is provided, make sure it's in array format (and stop single names becoming zero-dim array)
    else:
        name_list = np.array(copy.deepcopy(name))
        if name_list.shape == ():
            name_list = np.array([name_list.tolist()])
    del(name)

    # Do final check that all input sequences are the right length
    if np.std([float(ra_list.size), float(dec_list.size), float(width_list.size), float(name_list.size)]) > 0:
        raise Exception('Input sequences of ra, dec, with, and name all need to be the same length')

    # If no outout directory specified, set to current working directory
    if out_dir == None:
        out_dir = os.getcwd()

    # Check that output directory exists
    if not os.path.exists(out_dir):
        raise Exception('Specified output directory does not exist')

    # Create temporary directory
    if temp_dir == None:
        temp_dir = os.path.join(out_dir,'Temp')

    # Check that temp directory exists, if it does, warn user that contents may be overwritten
    if os.path.exists(temp_dir):
        print('Specificed temporary directory already exists; note that any existing contents may be overwritten')

    # Else, if temp directory doesn't already exist, create it
    else:
        os.mkdir(temp_dir)

    # Define dictionary of band properties
    bands_dict = {'Planck030':{'band_name':'Planck 030 I','freq':'30GHz','wavelength':'10600','pix_size':204.0,'beam_area':1188.945,'units':'K(cmb)','instrument':'LFI','conversion':0.979328*24.845597},
                  'Planck044':{'band_name':'Planck 044 I','freq':'44GHz','wavelength':'6810','pix_size':204.0,'beam_area':832.168,'units':'K(cmb)','instrument':'LFI','conversion':0.95121302*59.666236},
                  'Planck070':{'band_name':'Planck 070 I','freq':'70GHz','wavelength':'4260','pix_size':204.0,'beam_area':200.519,'units':'K(cmb)','instrument':'LFI','conversion':0.88140690*151.73238},
                  'Planck100':{'band_name':'Planck 100 I','freq':'100GHz','wavelength':'3000','pix_size':204.0,'beam_area':105.777,'units':'K(cmb)','instrument':'HFI','conversion':0.76581996*306.81118},
                  'Planck143':{'band_name':'Planck 143 I','freq':'143GHz','wavelength':'2100','pix_size':102.0,'beam_area':59.952,'units':'K(cmb)','instrument':'HFI','conversion':0.59714682*627.39818},
                  'Planck217':{'band_name':'Planck 217 I','freq':'217GHz','wavelength':'1380','pix_size':102.0,'beam_area':28.426,'units':'K(cmb)','instrument':'HFI','conversion':0.31573332*1444.7432},
                  'Planck353':{'band_name':'Planck 353 I','freq':'353GHz','wavelength':'850','pix_size':102.0,'beam_area':26.653,'units':'K(cmb)','instrument':'HFI','conversion':0.071041398*3823.1434},
                  'Planck545':{'band_name':'Planck 545 I','freq':'545GHz','wavelength':'550','pix_size':102.0,'beam_area':26.305,'units':'MJy/sr','instrument':'HFI'},
                  'Planck857':{'band_name':'Planck 857 I','freq':'857GHz','wavelength':'350','pix_size':102.0,'beam_area':23.985,'units':'MJy/sr','instrument':'HFI'}}

    # Record time taken
    time_list = [time.time()]

    # Loop over each source
    for i in np.random.permutation(range(name_list.shape[0])):
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        width = width_list[i]

        # If we're not repeating already-processed targets, check if this target has already been completed
        if not replace:
            bands_done = 0
            for band in bands_dict.keys():
                if os.path.exists(os.path.join(out_dir,name+'_Planck_'+bands_dict[band]['wavelength']+'.fits')):
                    bands_done += 1

                # Also check for null files, indicated data not available for a givne band
                elif os.path.exists(os.path.join(out_dir,'.'+name+'_Planck_'+bands_dict[band]['wavelength']+'.null')):
                    bands_done += 1

            # If this source has already been processed in all bands, skip it
            if bands_done == len(bands_dict.keys()):
                print('Planck data for '+name+ ' already generated; continuing to next target')
                time_list.append(time.time())
                continue
        print('Processing Planck data for target '+name)

        # In parallel, retrieve Planck data in each band from NASA SkyView
        """pool = mp.Pool(processes=5)"""
        for band in bands_dict.keys():
            #pool.apply_async( Planck_SkyView, args=(name, ra, dec, width, band, bands_dict, temp_dir,) )
            Planck_SkyView(name, ra, dec, width, band, bands_dict, temp_dir)
        """pool.close()
        pool.join()"""

        # In parallel, generate final standardised maps for each band
        """pool = mp.Pool(processes=9)"""
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( Planck_Generator, args=(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails,) )
            Planck_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails)
        """pool.close()
        pool.join()"""

        # Clean memory, and return timings (if more than one target being processed)
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        if len(name) > 1:
            print('Estimated time until Planck data completed for all targets: '+time_est)

    # Report completion
    print('Total time elapsed: '+str((time.time()-time_list[0])/3600.0)+' hours')

    # Tidy up (best as we can)
    gc.collect()
    try:
        shutil.rmtree(temp_dir)
    except:
        ChrisFuncs.RemoveCrawl(temp_dir)
        print('Unable to fully tidy up temporary directory; probably due to NFS locks on network drive')





# Define function to query for, and retrieve, Planck data from NASA SkyView
def Planck_SkyView(name, ra, dec, width, band, bands_dict, temp_dir):
    print('Querying for Planck '+bands_dict[band]['wavelength']+'um data for '+name+' from NASA SkyView')
    position_string = str(ra)+' '+str(dec)
    query_success = None
    pix_size = bands_dict[band]['pix_size']

    # Add tiny permutation to width to stop SkyView flailing over repeat queries
    width *= 1.0 + (0.01 * (np.random.rand() - 0.5))

    # Perform query
    while query_success==None:
        try:
            query_filename = os.path.join(temp_dir, name+'_Planck_'+bands_dict[band]['wavelength']+'.fits')
            query_url = SkyView.get_image_list(position_string, bands_dict[band]['band_name'], deedger='_skip_', pixels=int((width*3600.0)/pix_size), radius=astropy.units.Quantity(width, unit='deg'))
            if len(query_url)!=1:
                pdb.set_trace()
            query_success = True
        except:
            query_success = False

    # Retrieve and verify data
    if query_success:
        print('Retrieving identified Planck '+bands_dict[band]['wavelength']+'um data for '+name)
        Planck_wget(str(query_url[0]), query_filename)
        try:
            astropy.io.fits.info(query_filename)
        except Exception:
            query_success = False
            pdb.set_trace()

    # If no data available, generate null file and report failure
    if not query_success:
        os.system('touch '+os.path.join(temp_dir,'.'+name+'_Planck_'+bands_dict[band]['wavelength']+'.null'))
        print('No '+band+' data for '+name)






# Define function to finalise Planck image of a given source in a given band
def Planck_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails):
    wavelength = band_dict['wavelength']
    print('Generating final standardised map of Planck '+wavelength+'um data for '+name)

    # If null file exists for this target in this band, copy it to final output directory
    if os.path.exists(os.path.join(temp_dir,'.'+name+'_Planck_'+band_dict['wavelength']+'.null')):
        shutil.copy(os.path.join(temp_dir,'.'+name+'_Planck_'+band_dict['wavelength']+'.null'),
                  os.path.join(out_dir,'.'+name+'_Planck_'+band_dict['wavelength']+'.null'))
    else:

        # Read in map
        in_img, in_hdr = astropy.io.fits.getdata(os.path.join(temp_dir,name+'_Planck_'+wavelength+'.fits'), header=True)
        in_wcs = astropy.wcs.WCS(in_hdr)
        pix_width_arcsec = 3600.0 * astropy.wcs.utils.proj_plane_pixel_scales(in_wcs).mean()
        out_img = in_img.copy()

        """# Set non-coverage pixels to be nan
        out_img[ np.where(out_img==0) ] = np.NaN"""

        # Calculate sr/pixel
        sr_per_sqarcsec = 2.3504E-11
        sr_per_pixels = sr_per_sqarcsec * pix_width_arcsec**2

        # Convert pixel units from K(cmb) to MJy/sr
        if band_dict['units']=='K(cmb)':
            out_img *= band_dict['conversion']

        # If desired, convert pixel units from MJy/sr to Jy/pix
        pix_unit = 'MJy/sr'
        if flux:
            out_img *= 1E6
            out_img *= sr_per_pixels
            pix_unit = 'Jy/pix'

        # Create standard header
        out_hdr = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        out_hdr.set('TARGET', name, 'Target source of this map')
        out_hdr.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        out_hdr.set('SIGUNIT', pix_unit, 'Unit of the map')
        out_hdr.set('TELESCOP', 'Planck', 'Telescope that made this observation')
        out_hdr.set('INSTMNT', band_dict['instrument'], 'Instrument used for this observation')
        out_hdr.set('FREQ', band_dict['freq'], 'Filter used for this observation')
        out_hdr.set('WVLNGTH', wavelength+'um', 'Wavelength of observation')
        out_hdr.set('MAPDATE', date, 'Date this cutout was made from the existing reduced data')
        out_hdr['SOFTWARE'] = 'This cutout was produced using the Ancillary Data Button, written by Chris Clark, available from' \
                           + ' https://github.com/Stargrazer82301/AncillaryDataButton/, following procedures laid out in' \
                           + ' Clark et al (2018, A&A 609 A37) and Saintonge et al (2018).'

        # Construct WCS system, and append to header
        cutout_wcs = astropy.wcs.WCS(naxis=2)
        cutout_wcs.wcs.crpix = [in_hdr['CRPIX1'], in_hdr['CRPIX2']]
        cutout_wcs.wcs.cdelt = [in_hdr['CDELT1'], in_hdr['CDELT2']]
        cutout_wcs.wcs.crval = [in_hdr['CRVAL1'], in_hdr['CRVAL2']]
        cutout_wcs.wcs.ctype = [in_hdr['CTYPE1'], in_hdr['CTYPE2']]
        cutout_wcs_header = cutout_wcs.to_header()
        for card in cutout_wcs_header.cards:
            out_hdr.append(card)

        # Write output FITS file
        astropy.io.fits.writeto(os.path.join(out_dir,name+'_Planck_'+wavelength+'.fits'), data=out_img, header=out_hdr, overwrite=True)

        # Make thumbnail image of cutout
        if thumbnails:
            try:
                thumb_out = aplpy.FITSFigure(out_dir+name+'_Planck_'+wavelength+'.fits')
                thumb_out.show_colorscale(cmap='gist_heat', stretch='arcsinh')
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=500, lw=2.5, edgecolor='#01DF3A')
                thumb_out.save(os.path.join(out_dir,name+'_Planck_'+wavelength+'.png'), dpi=125)
                thumb_out.close()
            except:
                print('Failed making thumbnail for '+name)
                pdb.set_trace()

        # Clean memory before finishing
        gc.collect()



# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget DSS tiles
def Planck_wget(tile_url, tile_filename):
    if os.path.exists(tile_filename):
        os.remove(tile_filename)
    success = False
    fail_count = 0
    while success==False:
        if fail_count == 4:
            return
        try:
            urllib.request.urlretrieve(tile_url, tile_filename)
            success = True
        except:
            fail_count += 1
            time.sleep(5.0)
            success = False



# Define function to provide appropriate padding for FITS header entires
def Padding(entry):
    return ( ' ' * np.max([ 0, 19-len( str(entry) ) ]) ) + ' / '