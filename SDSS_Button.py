# Import smorgasbord
import pdb
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import multiprocessing as mp
import gc
import re
import copy
import shutil
import ssl
import time
import datetime
import signal
import astropy.io.fits
import astropy.wcs
import astropy.io.votable
import aplpy
import montage_wrapper
import wget
import ChrisFuncs
plt.ioff()





# Define main function
def Run(ra, dec, width, name=None, out_dir=None, temp_dir=None, replace=False, flux=True, thumbnails=False):
    """
    Function to generate standardised cutouts of SDSS observations.

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
                A string giving the path to be used as a temporary working directory by SDSS_Button. If not provided,
                a temporary directory will be created inside the output directory.
        replace: bool, optional
                If False, SDSS_Button will search the output directory for any pre-existing output FITS files from
                previous runs of the function, and will not bother repeat creating these maps (making it easy to resume
                processing a large number of targets from an interruption. If True, SDSS_Button will produce maps for
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

    # Read in list of SDSS DR12 primary fields
    dr12_pri = np.genfromtxt('SDSS_DR12_Primary_Fields.dat.gz', names=True)
    dr12_pri = [ str(int(pri['RUN']))+' '+str(int(pri['CAMCOL']))+' '+str(int(pri['FIELD'])) for pri in dr12_pri ]

    # Get around the fact that the SDSS FTP server doesn't have valid SSl certificates
    ssl._create_default_https_context = ssl._create_unverified_context

    # Define dictionary of band properties
    bands_dict = {'u':{'band_name':'u','wavelength':'255.10','sdss_to_ab':-0.04,},
                  'g':{'band_name':'g','wavelength':'468.60','sdss_to_ab':0.0,},
                  'r':{'band_name':'r','wavelength':'616.60','sdss_to_ab':0.0,},
                  'i':{'band_name':'i','wavelength':'748.00','sdss_to_ab':0.0,},
                  'z':{'band_name':'z','wavelength':'893.00','sdss_to_ab':+0.02,}}

    # Record time taken
    time_list = [time.time()]

    # Loop over each target
    for i in np.random.permutation(range(name_list.shape[0])):
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        width = width_list[i]

        # If we're not repeating already-processed targets, check if this target has already been completed
        if not replace:
            bands_done = 0
            for band in bands_dict.keys():
                if os.path.exists(os.path.join(out_dir,name+'_SDSS_'+bands_dict[band]['wavelength']+'.fits')):
                    bands_done += 1

                # Also check for null files, indicated data not available for a givne band
                elif os.path.exists(os.path.join(out_dir,'.'+name+'_SDSS_'+bands_dict[band]['wavelength']+'.null')):
                    bands_done += 1

            # If this source has already been processed in all bands, skip it
            if bands_done == len(bands_dict.keys()):
                print('SDSS data for '+name+ ' already processed (if available); continuing to next target')
                time_list.append(time_list)
                continue
        print('Processing SDSS data for target '+name)

        # Create field processing dirctories (deleting any prior), and set appropriate Python (ie, SWarp) working directory
        gal_dir = os.path.join(temp_dir,str(name))+'/'
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)
        os.makedirs(os.path.join(gal_dir,'Raw'))
        os.chdir(os.path.join(gal_dir,'Raw'))

        # Use Montage to perform query of SDSS DR9 fields (inluding additional 0.14 degrees to search box to account for field size?)
        coverage = True
        bands_null = 0
        sdss_urls = []
        for band in bands_dict.keys():
            if coverage==False:
                continue
            montage_wrapper.commands.mArchiveList('SDSS', band.lower(), str(ra)+' '+str(dec), width, width, gal_dir+'Query_'+band+'.txt')#width*(2.**0.5), width*(2.**0.5)

            # Check if query returned any results
            if os.stat(gal_dir+'Query_'+band+'.txt').st_size==0:
                bands_null += 1
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_SDSS_'+band+'.null'))
                continue
            elif sum(1 for line in open(gal_dir+'Query_'+band+'.txt'))<=3:
                bands_null += 1
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_SDSS_'+band+'.null'))
                continue

            # If SDSS coverage present, extract field info from query results, removing duplicate entries
            else:
                montage_urls = np.genfromtxt(gal_dir+'Query_'+band+'.txt', skip_header=2, usecols=[9], dtype=[('URL','S256')])
                [ sdss_urls.append(montage_url[0]) for montage_url in montage_urls ]



        # Remove duplicate URLs from URL list, and only use primary runs
        sdss_urls = list(set(sdss_urls))
        sdss_urls_pri = SDSS_Primary_Check(sdss_urls, dr12_pri)
        print(str(len(sdss_urls_pri))+' of '+str(len(sdss_urls))+' urls for '+name+' are primary')
        sdss_urls = sdss_urls_pri
        if len(sdss_urls)==0:
            bands_null = len(bands_dict.keys())

        # Update URLS to DR12 location
        sdss_urls = [sdss_url.replace('dr10','dr12') for sdss_url in sdss_urls]

        # If no SDSS coverge, and skip onwards
        if (coverage == False) or (bands_null == len(bands_dict.keys())):
            print('No SDSS coverage for '+name+'; continuing to next target')
            time_list.append(time_list)
            continue

        # In parallel, download SDSS fields (NB: SDSS server will not permit more than 5 simultaneous wget downloads)
        print('Downloading '+str(len(sdss_urls))+' SDSS fields for '+name)
        dl_complete = False
        dl_fail_count = 0
        while not dl_complete:
            try:
                signal.alarm(7200)
                dl_pool = mp.Pool(processes=5)
                for sdss_url in sdss_urls:
                    dl_pool.apply_async( SDSS_wget, args=(sdss_url, gal_dir+'Raw/',) )
                    #SDSS_wget(sdss_url, gal_dir+'Raw/')
                dl_pool.close()
                dl_pool.join()
                dl_complete = True
                signal.alarm(0)
            except:
                dl_fail_count += 1
                gc.collect()
                shutil.rmtree(os.path.join(gal_dir,'Raw'))
                os.makedirs(os.path.join(gal_dir,'Raw'))
                os.chdir(os.path.join(gal_dir,'Raw'))
                if dl_fail_count==5:
                    dl_complete = True
                print('Download sequence failed; reattemping')

        # In parallel, extract SDSS fields from their bz2 archives
        raw_files = np.array(os.listdir(gal_dir+'Raw/'))
        bz2_pool = mp.Pool(processes=mp.cpu_count())
        for raw_file in raw_files:
            bz2_pool.apply_async( SDSS_Extract, args=(gal_dir+'Raw/'+raw_file,) )
        bz2_pool.close()
        bz2_pool.join()

        # Perform metadata query to find if any fields provide coverage at particular coord of target; if not, clean up and move on
        montage_wrapper.commands.mImgtbl(os.path.join(gal_dir,'Raw'),  gal_dir+'Metadata_'+band+'.dat', corners=True)
        montage_wrapper.commands_extra.mCoverageCheck(gal_dir+'Metadata_'+band+'.dat', gal_dir+'Metada_Query_'+band+'.dat', ra=ra, dec=dec, mode='point')
        if os.stat(gal_dir+'Metada_Query_'+band+'.dat').st_size>0:
            pass
        elif sum(1 for line in open(gal_dir+'Metada_Query_'+band+'.dat'))>=3:
            pass
        else:
            print(name+' not covered by SDSS')
            shutil.rmtree(gal_dir)
            continue

        # Copy files to relevant folders
        for band in bands_dict.keys():
            if os.path.exists(gal_dir+band+'/'):
                shutil.rmtree(gal_dir+band+'/')
            os.makedirs(gal_dir+band+'/')
            raw_files = np.array(os.listdir(gal_dir+'Raw/'))
            for raw_file in raw_files:
                if '-'+band+'-' in raw_file:
                    shutil.move(gal_dir+'Raw/'+raw_file, gal_dir+band)

        # In parallel, Montage together each band's input fields, whilst dealing with timeouts
        out_pix_width = 0.396
        complete = False
        while not complete:
            try:
                #signal.alarm(3600)
                pool = mp.Pool(processes=5)
                for band in bands_dict.keys():
                    input_dir = gal_dir+band+'/'
                    pool.apply_async( SDSS_Montage, args=(name, ra, dec, out_pix_width, width, band, input_dir, temp_dir,) )
                    #SDSS_Montage(name, ra, dec, out_pix_width, width, band, input_dir, temp_dir)
                pool.close()
                pool.join()
                complete = True
                signal.alarm(0)
            except:
                gc.collect()
                if os.path.exists(os.path.join(input_dir,'Montage_Temp')):
                    shutil.rmtree(os.path.join(input_dir,'Montage_Temp'))
                print('Mosaicing failed')

        # In parallel, generate final standardised maps for each band
        pool = mp.Pool(processes=9)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( SDSS_Generator, args=(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails,) )
            SDSS_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails)
        pool.close()
        pool.join()

        # Clean memory, and return timings (if more than one target being processed)
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        if len(name) > 1:
            print('Estimated time until SDSS data completed for all targets: '+time_est)

    # Report completion
    try:
        shutil.rmtree(temp_dir)
    except:
        try:
            time.sleep(10.0)
            shutil.rmtree(temp_dir)
        except:
            pdb.set_trace()
    gc.collect()
    print('All available SDSS imagery acquired for all targets')





# Define function to Montage together contents of folder
def SDSS_Montage(name, ra, dec, pix_width, map_width, band, input_dir, out_dir):
    print('Montaging '+name+'_SDSS_'+band)
    location_string = str(ra)+' '+str(dec)
    if os.path.exists(input_dir+'Montage_Temp'):
        shutil.rmtree(input_dir+'Montage_Temp')
    os.makedirs(input_dir+'Montage_Temp')
    os.chdir(input_dir+'Montage_Temp')
    montage_wrapper.commands.mHdr(location_string, map_width, input_dir+'Montage_Temp/'+str(name)+'_HDR', pix_size=pix_width)
    montage_wrapper.commands.mExec('SDSS', band.lower(), raw_dir=input_dir, level_only=False, debug_level=0, output_image=os.path.join(out_dir,name+'_SDSS_'+band+'.fits'), region_header=input_dir+'Montage_Temp/'+str(name)+'_HDR', workspace_dir=input_dir+'Montage_Temp')
    shutil.rmtree(input_dir+'Montage_Temp')
    print('Completed Montaging '+name+'_SDSS_'+band)



# Define function to finalise SDSS image of a given source in a given band
def SDSS_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails):
    band = band_dict['band_name']
    wavelength = band_dict['wavelength']
    print('Generating final standardised map of SDSS '+band+' data for '+name)

    # If null file exists for this target in this band, copy it to final output directory
    if os.path.exists(os.path.join(temp_dir,'.'+name+'_SDSS_'+band+'.null')):
        shutil.copy(os.path.join(temp_dir,'.'+name+'_SDSS_'+band+'.null'),
                  os.path.join(out_dir,'.'+name+'_SDSS_'+band+'.null'))
    else:

        # Read in map
        in_img, in_hdr = astropy.io.fits.getdata(os.path.join(temp_dir,name+'_SDSS_'+band+'.fits'), header=True)
        in_wcs = astropy.wcs.WCS(in_hdr)
        in_pix_width_arcsec = 3600.0 * np.abs( np.max( in_wcs.pixel_scale_matrix ) )
        out_img = in_img.copy()

        # Add minimum pixel value to image, to prevent NaNs from appearing later
        out_img = out_img + np.abs(2.0 * np.nanmin(out_img))

        # Convert each pixel value to SDSS magnitudes
        out_img = 22.5 - (2.5 * np.log10(out_img))

        # Convert each pixel value to AB magnitudes
        out_img = out_img + band_dict['sdss_to_ab']

        # Convert each pixel value to janskys
        out_img = ChrisFuncs.ABMagsToJy(out_img)

        # Scale pixel values to account pixel size increase from 0.396 to 0.45 arcseconds (as Montage conserves surface brightness)
        out_img *= in_pix_width_arcsec**2 / 0.396**2

        # If desired, convert pixel units from Jy/pix to Jy/pix
        pix_unit = 'Jy/pix'
        if not flux:
            sr_per_sqarcsec = 2.3504E-11
            sr_per_pixels = sr_per_sqarcsec * in_pix_width_arcsec**2
            out_img *= 1E6
            out_img *= sr_per_pixels
            pix_unit = 'MJy/sr'

        # Create standard header
        out_hdr = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        out_hdr.set('TARGET', name, 'Target source of this map')
        out_hdr.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        out_hdr.set('SIGUNIT', pix_unit, 'Unit of the map')
        out_hdr.set('TELESCOP', 'SDSS', 'Telescope that made this observation')
        out_hdr.set('FILTER', band, 'Filter used for this observation')
        out_hdr.set('WVLNGTH', wavelength+'nm', 'Effective wavelength of observation')
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
        astropy.io.fits.writeto(os.path.join(out_dir,name+'_SDSS_'+band+'.fits'), data=out_img, header=out_hdr, overwrite=True)

        # Make thumbnail image of cutout
        if thumbnails:
            try:
                thumb_out = aplpy.FITSFigure(out_dir+name+'_SDSS_'+band+'.fits')
                thumb_out.show_colorscale(cmap='gist_heat', stretch='arcsinh')
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=500, lw=2.5, edgecolor='#01DF3A')
                thumb_out.save(os.path.join(out_dir,name+'_SDSS_'+band+'.jpg'), dpi=125)
                thumb_out.close()
            except:
                print('Failed making thumbnail for '+name)
                pdb.set_trace()

        # Clean memory before finishing
        gc.collect()


# Define function to wget SDSS fields; reject and re-acquire fields less than 1MB in size
def SDSS_wget(frame_url, path):
    fitsname = frame_url.split('/')[-1:][0]
    print('Acquiring '+fitsname)
    if os.path.exists(path+fitsname):
        os.remove(path+fitsname)
    success = False
#    while success==False:
#        try:
    wget.download(frame_url, out=path+fitsname)
    print(path+fitsname)
    filesize = os.stat(path+fitsname).st_size
    if float(filesize)<1048576.0:
        raise NameError('File not large enough')
    print('Successful acquisition of '+fitsname)
    success = True
#        except:
#            print('Failure! Retrying acquistion of '+fitsname)
#            time.sleep(0.1)
#            success = False



# Define function to extract SDSS bz2 archives
def SDSS_Extract(filepath):
    #print('Decompressing file '+str( filepath.split('/')[-1:][0] ))
    os.system('bzip2 -d '+filepath)



# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to check if a list of SDSS DR12 URL correspond to a primary fields, and return only those that are
def SDSS_Primary_Check(urls, index):
    urls_pri = []
    for url in urls:
        url = url.decode('utf-8')
        run = url.split('/')[9].lstrip('0')
        camcol = url.split('/')[10].lstrip('0')
        field = url.split('/')[11].split('.')[0].split('-')[4].lstrip('0')
        check_string = run+' '+camcol+' '+field
        if check_string in index:
            urls_pri.append(url)
    return urls_pri



# Define function to consturct URL for a given SDSS field
def SDSS_URL(run, camcol, field, band):
    sdss_url = 'http://data.sdss3.org/sas/dr12/boss/photoObj/frames/301/'+run+'/'+camcol+'/frame-'+band+'-'+run.zfill(6)+'-'+camcol+'-'+field.zfill(4)+'.fits.bz2'
    return sdss_url




Run(023.462100, 30.659942, 2.0, name='M33', out_dir='/astro/dust_kg/cclark/Lea/')