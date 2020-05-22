# Import smorgasbord
import pdb
import os
import sys
import multiprocessing as mp
import numpy as np
import astropy.io.votable
import astropy.coordinates
import reproject
import warnings
warnings.filterwarnings("ignore")
import aplpy
import urllib
import gc
import re
import shutil
import copy
import time
import datetime
import joblib
from ChrisFuncs import TimeEst, RemoveCrawl
from ChrisFuncs.Fits import FitsHeader



# Define main function
def Run(ra, dec, width, name=None, out_dir=None, temp_dir=None, replace=False, flux=True, thumbnails=False, montage_path=None):
    """
    Function to generate standardised cutouts of IRAS-ISSA observations from the calibrated plates hosted on IRSA.

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
                A string giving the path to be used as a temporary working directory by ISSA_Button. If not provided,
                a temporary directory will be created inside the output directory.
        replace: bool, optional
                If False, ISSA_Button will search the output directory for any pre-existing output FITS files from
                previous runs of the function, and will not bother repeat creating these maps (making it easy to resume
                processing a large number of targets from an interruption. If True, ISSA_Button will produce maps for
                all input targets, regardless of whether maps for these targets already exist in the output directory.
        flux: bool, optional
                If True, output maps will be in flux density units of Jy/pix. If false, output maps will be in surface
                brightness units of MJy/sr.
        thumbnails: bool, optional
                If True, JPG thumbnail images of the generated maps will also be proced and placed in out_dir.
        montage_path
                This function requires Montage to be installed. If your Montage commands are found in a location that
                the montage_wapper package will not find by default, provide a string here that gives the path to the
                directory the commands are found in.
    """

    # If Montage commands directory provided, append it to path
    try:
        if montage_path != None:
            sys.path.append(montage_path)
            os.environ['PATH'] = os.environ['PATH'] + ':' + montage_path
        import montage_wrapper
    except:
        if montage_path != None:
            sys.path.append(montage_path)
            os.environ['PATH'] = os.environ['PATH'] + ':' + montage_path
        import montage_wrapper
    montage_wrapper = montage_wrapper

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
    bands_dict = {'12':{'band_name':'ISSA  12','wavelength':'12','band_num':'1','pix_size':90.0},
                  '25':{'band_name':'ISSA  25','wavelength':'25','band_num':'2','pix_size':90.0},
                  '60':{'band_name':'ISSA  60','wavelength':'60','band_num':'3','pix_size':90.0},
                  '100':{'band_name':'ISSA 100','wavelength':'100','band_num':'4','pix_size':90.0}}

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
                if os.path.exists(os.path.join(out_dir,name+'_ISSA_'+bands_dict[band]['wavelength']+'.fits')):
                    bands_done += 1

                # Also check for null files, indicated data not available for a givne band
                elif os.path.exists(os.path.join(out_dir,'.'+name+'_ISSA_'+bands_dict[band]['wavelength']+'.null')):
                    bands_done += 1

            # If this source has already been processed in all bands, skip it
            if bands_done == len(bands_dict.keys()):
                print('IRAS-ISSA data for '+name+ ' already generated; continuing to next target')
                time_list.append(time.time())
                continue
        print('Processing IRAS-ISSA data for target '+name)

        # In parallel, retrieve IRAS-ISSA data in each band from IRSA
        pool = mp.Pool(processes=4)
        for band in bands_dict.keys():
            pool.apply_async(ISSA_Query, args=(name, ra, dec, width, band, bands_dict, temp_dir, montage_path,))
            #ISSA_Query(name, ra, dec, width, band, bands_dict, temp_dir, montage_path=montage_path)
        pool.close()
        pool.join()

        # Generate final standardised maps for each band (this can be run in parallel, but that seems to cause thumbnail problems)
        pool = mp.Pool(processes=4)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async(ISSA_Generator, args=(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails,))
            ISSA_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails)
        pool.close()
        pool.join()

        # Clean memory, and return timings (if more than one target being processed)
        gc.collect()
        time_list.append(time.time())
        time_est = TimeEst(time_list, len(name_list))
        if len(name) > 1:
            print('Estimated time until IRAS-ISSA data completed for all targets: '+time_est)

        # Tidy up (best as we can)
        gc.collect()
        try:
            shutil.rmtree(temp_dir)
        except:
            RemoveCrawl(temp_dir)
            print('Unable to fully tidy up temporary directory; probably due to NFS locks on network drive')

    # Report completion
    print('Total time elapsed: '+str((time.time()-time_list[0])/3600.0)+' hours')



# Define function to retrieve, select, and mosaic IRAS-ISSA data from IRSA
def ISSA_Query(name, ra, dec, width, band, bands_dict, temp_dir, montage_path=None):

    # If Montage commands directory provided, append it to path
    try:
        import montage_wrapper
    except:
        sys.path.append(montage_path)
        os.environ['PATH'] = os.environ['PATH'] + ':' + montage_path
        import montage_wrapper

    # Generate list of all ISSA plate fields in this band (which take form iYYYBXh0.fits, where YYY is a number between 001 and 430, and X is the field between 1 and 4)
    issa_url = 'https://irsa.ipac.caltech.edu/data/ISSA/ISSA_complete_v2/'
    issa_fields = np.arange(1,431).astype(str)
    issa_fields = [''.join(['i',field.zfill(3),'bXh0']) for field in issa_fields]

    # Check if a folder for the raw ISSA plates exists in the temporary directory; if not, create it
    print('Ensuring all raw '+bands_dict[band]['wavelength']+'um IRAS-ISSA plates are available')
    band = bands_dict[band]['wavelength']
    raw_dir = os.path.join(temp_dir,'Raw',band)
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)

    # Look to see if all ISSA fields for this band are already present in the temporary directory; if not, wget them
    wget_list = []
    for issa_field in np.random.permutation(issa_fields):
        issa_ref_file = issa_field.replace('X',bands_dict[band]['band_num'])+'.fit'
        issa_ref_path = os.path.join(raw_dir,issa_ref_file)
        if not os.path.exists(issa_ref_path):
            wget_list.append([issa_url+issa_ref_file,issa_ref_path])
    if len(wget_list) > 0:
        print('Downloading raw '+bands_dict[band]['wavelength']+'um IRAS-ISSA plates (note that this will entail downloding up to ~4GB of data)')
        if mp.current_process().name == 'MainProcess':
            joblib.Parallel( n_jobs=mp.cpu_count()-1 )\
                           ( joblib.delayed( ISSA_wget )\
                           ( wget_list[w][0], wget_list[w][1] )\
                           for w in range(len(wget_list)) )
        else:
            for w in range(len(wget_list)):
                os.system('curl '+wget_list[w][0]+' -o '+'"'+wget_list[w][1]+'"')


    # If image metadata table doesn't yet exist for this band, run mImgtbl over raw data to generate it
    mImgtbl_tablepath = os.path.join(raw_dir,'ISSA_'+band+'_Metadata_Table.tbl')
    if not os.path.exists(mImgtbl_tablepath):
        montage_wrapper.mImgtbl(raw_dir, mImgtbl_tablepath, corners=True)

    # Now that we know we have data, set up processing for this source in particular
    print('Computing overlap of '+bands_dict[band]['wavelength']+'um IRAS-ISSA plates with '+name)
    ra, dec, width = float(ra), float(dec), float(width)
    pix_size = bands_dict[band]['pix_size']

    # Find which plates have coverage over our target region
    mCoverageCheck_tablepath = os.path.join(raw_dir,u'ISSA_'+band+'_Coverage_Table.tbl')
    if os.path.exists(mCoverageCheck_tablepath):
        os.remove(mCoverageCheck_tablepath)
    montage_wrapper.mCoverageCheck(mImgtbl_tablepath, mCoverageCheck_tablepath, ra=ra, dec=dec, mode='box', width=width)

    # Read in coveage tables to identify what plates we need; if no coverage, write null output file and stop here
    print('Reprojecting IRAS-ISSAS '+bands_dict[band]['wavelength']+'um plates that cover '+name)
    mCoverageCheck_table = np.genfromtxt(mCoverageCheck_tablepath, skip_header=3, dtype=None, encoding=None)
    if mCoverageCheck_table.size == 0:
        os.system('touch '+os.path.join(temp_dir,'.'+name+'_IRAS-IRIS_'+band+'.null'))
        print('No IRAS-IRIS '+band+'um data for '+name)
        return
    reproj_dir = os.path.join(temp_dir,'Reproject',band)
    if not os.path.exists(reproj_dir):
        os.makedirs(reproj_dir)

    # Extract paths from coverage table, with handling for weird astropy behavior when table has only one row
    if mCoverageCheck_table.size == 1:
        raw_paths = [mCoverageCheck_table['f36'].tolist()]
    else:
        raw_paths = [str(mCoverageCheck_table['f36'][i]) for i in range(mCoverageCheck_table['f36'].size)]
    reproj_paths = [raw_paths[i].replace(raw_dir,reproj_dir) for i in range(len(raw_paths))]
    reproj_hdr = FitsHeader(ra, dec, width, pix_size)

    # Reproject identified plates in turn (dealing with possible corrupt downloads, and stupid unecessary third axis, grrr)
    for i in range(len(raw_paths)):
        raw_path, reproj_path = raw_paths[i], reproj_paths[i]
        try:
            raw_img, raw_hdr = astropy.io.fits.getdata(raw_path, header=True, memmap=False)
        except:
            raw_url = issa_url + raw_path.split('/')[-1]
            os.system('curl '+raw_url+' -o '+'"'+raw_path+'"')
            raw_img, raw_hdr = astropy.io.fits.getdata(raw_path, header=True, memmap=False)
        raw_hdr.set('NAXIS', value=2)
        raw_hdr.remove('NAXIS3')
        raw_hdr.remove('CRVAL3')
        raw_hdr.remove('CRPIX3')
        raw_hdr.remove('CTYPE3')
        raw_hdr.remove('CDELT3')
        raw_img = raw_img[0,:,:]
        raw_hdu = astropy.io.fits.PrimaryHDU(data=raw_img, header=raw_hdr)
        reproj_img = reproject.reproject_exact(raw_hdu, reproj_hdr, parallel=False)[0]
        astropy.io.fits.writeto(reproj_path, data=reproj_img, header=reproj_hdr, overwrite=True)
        del(raw_hdu)
        del(raw_img)
        del(raw_hdr)
        gc.collect()

    # Now mosaic the reprojected images
    mosaic_list = []
    [mosaic_list.append(astropy.io.fits.getdata(reproj_path)) for reproj_path in reproj_paths]
    mosaic_array = np.array(mosaic_list)
    mosaic_img = np.nanmean(mosaic_array, axis=0)
    mosaic_hdr = FitsHeader(ra, dec, width, pix_size)

    """# Write finished mosaic to file
    astropy.io.fits.writeto(os.path.join(temp_dir,name+'_IRAS-ISSA_'+band+'.fits'), data=mosaic_img, header=mosaic_hdr, overwrite=True)"""

    # Check that target coords have coverage in mosaic
    mosaic_wcs = astropy.wcs.WCS(mosaic_hdr)
    mosaic_centre = mosaic_wcs.all_world2pix([[ra]], [[dec]], 0, ra_dec_order=True)
    mosaic_i, mosaic_j = mosaic_centre[1][0], mosaic_centre[0][0]
    if np.isnan(mosaic_img[int(np.round(mosaic_i)),int(np.round(mosaic_j))]):
        os.system('touch '+os.path.join(temp_dir,'.'+name+'_IRAS-ISSA_'+band+'.null'))
        print('No IRAS-ISSA '+band+'um data for '+name)

    # If mosaic is good, write it to temporary directory
    else:
        astropy.io.fits.writeto(os.path.join(temp_dir,name+'_IRAS-ISSA_'+band+'.fits'), data=mosaic_img, header=mosaic_hdr, overwrite=True)



# Define function to finalise IRAS-ISSA image of a given source in a given band
def ISSA_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails):
    wavelength = band_dict['wavelength']
    print('Generating final standardised map of IRAS-ISSA '+wavelength+'um data for '+name)

    # If null file exists for this target in this band, copy it to final output directory
    if os.path.exists(os.path.join(temp_dir,'.'+name+'_IRAS-ISSA_'+band_dict['wavelength']+'.null')):
        shutil.copy(os.path.join(temp_dir,'.'+name+'_IRAS-ISSA_'+band_dict['wavelength']+'.null'),
                  os.path.join(out_dir,'.'+name+'_IRAS-ISSA_'+band_dict['wavelength']+'.null'))
    else:

        # Read in map
        in_img, in_hdr = astropy.io.fits.getdata(os.path.join(temp_dir,name+'_IRAS-ISSA_'+wavelength+'.fits'), header=True)
        in_wcs = astropy.wcs.WCS(in_hdr)
        pix_width_arcsec = 3600.0 * astropy.wcs.utils.proj_plane_pixel_scales(in_wcs).mean()
        out_img = in_img.copy()

        # Calculate sr/pixel
        sr_per_sqarcsec = 2.3504E-11
        sr_per_pixels = sr_per_sqarcsec * pix_width_arcsec**2

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
        out_hdr.set('TELESCOP', 'IRAS', 'Telescope that made this observation')
        out_hdr.set('INSTRMNT', 'IRAS', 'Instrument used for this observation')
        out_hdr.set('PIPELINE', 'ISSA (IRAS Sky Survey Atlas)', 'Data products from which this cutout was produced')
        out_hdr.set('WVLNGTH', wavelength+'um', 'Wavelength of observation')
        out_hdr.set('MAPDATE', date, 'Date this cutout was made from the existing reduced data')
        out_hdr.set('SOFTWARE', 'The Ancillary Data Button',
                    'This cutout was produced using the Ancillary Data Button, written by Chris Clark, available from' \
                    + ' https://github.com/Stargrazer82301/AncillaryDataButton/, following procedures laid out in' \
                    + ' Clark et al (2018, A&A 609 A37) and Saintonge et al (2018).')

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
        astropy.io.fits.writeto(os.path.join(out_dir,name+'_IRAS-ISSA_'+wavelength+'.fits'), data=out_img, header=out_hdr, overwrite=True)

        # Make thumbnail image of cutout
        if thumbnails:
            print('Generating thumbnail image of IRAS-ISSA '+wavelength+'um data for '+name)
            thumb_out = aplpy.FITSFigure(out_dir+name+'_IRAS-ISSA_'+wavelength+'.fits')
            thumb_out.show_colorscale(cmap='gist_heat', stretch='arcsinh')
            thumb_out.axis_labels.hide()
            thumb_out.tick_labels.hide()
            thumb_out.ticks.hide()
            thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=500, lw=2.5, edgecolor='#01DF3A')
            thumb_out.save(os.path.join(out_dir,name+'_IRAS-ISSA_'+wavelength+'.png'), dpi=125)
            thumb_out.close()

        # Clean memory before finishing
        gc.collect()



# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget DSS tiles
def ISSA_wget(tile_url, tile_filename):
    if os.path.exists(tile_filename):
        os.remove(tile_filename)
    success = False
    while success==False:
        try:
            urllib.request.urlretrieve(tile_url, tile_filename)
            success = True
        except:
            time.sleep(0.1)
            success = False


# Define function to provide appropriate padding for FITS header entires
def Padding(entry):
    return ( ' ' * np.max([ 0, 19-len( str(entry) ) ]) ) + ' / '