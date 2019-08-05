# Import smorgasbord
import pdb
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import multiprocessing as mp
import gc
import re
import copy
import shutil
import lmfit
import joblib
import time
import datetime
import signal
import astropy.io.fits
import astropy.wcs
import astropy.io.votable
import reproject
import aplpy
import wget
import ChrisFuncs
import ChrisFuncs.Fits
import ChrisFuncs.Coadd
plt.ioff()





# Define main function
def Run(ra, dec, width, name=None, out_dir=None, temp_dir=None, replace=False, flux=True, thumbnails=False, swarp_bg=False):
    """
    Function to generate standardised cutouts of GALEX observations.

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
                A string giving the path to be used as a temporary working directory by galex_Button. If not provided,
                a temporary directory will be created inside the output directory.
        replace: bool, optional
                If False, galex_Button will search the output directory for any pre-existing output FITS files from
                previous runs of the function, and will not bother repeat creating these maps (making it easy to resume
                processing a large number of targets from an interruption. If True, galex_Button will produce maps for
                all input targets, regardless of whether maps for these targets already exist in the output directory.
        flux: bool, optional
                If True, output maps will be in flux density units of Jy/pix. If false, output maps will be in surface
                brightness units of MJy/sr.
        thumbnails: bool, optional
                If True, JPG thumbnail images of the generated maps will also be proced and placed in out_dir.
        swarp_bg: bool, optional
                As GALEX_Button uses SWarp to perform the coaddition of the data, this kwarg provides the option to
                enable SWarp's background-sutraction (disabled be default, as it can remove large-scale structure).
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

    bands_dict = {'FUV':{'band_short':'fd','band_long':'FUV','wavelength':'152.8','zero_point':18.82},
                  'NUV':{'band_short':'nd','band_long':'NUV','wavelength':'227.1','zero_point':20.08}}

    # Provide information about different map types
    map_suffixes = ['-int.fits.gz','-skybg.fits.gz','-exp.fits.gz','-rr.fits.gz']

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
                if os.path.exists(os.path.join(out_dir,name+'_GALEX_'+bands_dict[band]['band_long']+'.fits')):
                    bands_done += 1

                # Also check for null files, indicated data not available for a givne band
                elif os.path.exists(os.path.join(out_dir,'.'+name+'_GALEX_'+bands_dict[band]['band_long']+'.null')):
                    bands_done += 1

            # If this source has already been processed in all bands, skip it
            if bands_done == len(bands_dict.keys()):
                print('GALEX data for '+name+ ' already processed (if available); continuing to next target')
                time_list.append(time.time())
                continue
        print('Processing GALEX data for target '+name)

        # Create field processing dirctories (deleting any prior), and set appropriate Python (ie, SWarp) working directory
        gal_dir = os.path.join(temp_dir,str(name))+'/'

        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)
        os.makedirs(os.path.join(gal_dir,'Download'))
        os.chdir(os.path.join(gal_dir,'Download'))

        # Perform query (removing pre-existing query file, if present)
        print('Querying GALEX MAST server using SIA protocol')
        query_url = 'https://mast.stsci.edu/portal/Mashup/VoQuery.asmx/SiaV1?MISSION=GALEX&POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS&FORMAT=image/fits'
        query_filename = os.path.join(temp_dir,name,str(name)+'.vot')
        if os.path.exists(query_filename):
            os.remove(query_filename)
        wget.download(query_url, out=query_filename)

        # Read query result VOTable
        query_table = astropy.table.Table.read(query_filename)

        # If query finds no matches, continue to next target
        if len(query_table)==0:
            print('No GALEX coverage for '+name+'; continuing to next target')
            time_list.append(time.time())
            for band in bands_dict.keys():
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_GALEX_'+band+'.null'))
            continue

        # Grab fits files URLs
        galex_urls = []
        for i in range(len(query_table)):
            galex_urls.append(query_table[i]['accessURL'])

        # In parallel, download GALEX fields
        print('Downloading '+str(len(galex_urls))+' GALEX fields for '+name)
        dl_complete = False
        dl_fail_count = 0
        while not dl_complete:
            try:
                signal.alarm(18000)
                dl_pool = mp.Pool(processes=20)
                for galex_url in galex_urls:
                    galex_url = str(galex_url, encoding='utf-8')
                    for map_suffix in map_suffixes:
                        map_url = galex_url.replace('-int.fits.gz', map_suffix)
                        dl_pool.apply_async( GALEX_wget, args=(map_url, os.path.join(gal_dir,'Download',map_url.split('/')[-1]),) )
                        #GALEX_wget(map_url, os.path.join(gal_dir,'Download',map_url.split('/')[-1]))
                dl_pool.close()
                dl_pool.join()
                dl_complete = True
                signal.alarm(0)
            except:
                dl_fail_count += 1
                gc.collect()
                shutil.rmtree(os.path.join(gal_dir,'Download'))
                os.makedirs(os.path.join(gal_dir,'Download'))
                os.chdir(os.path.join(gal_dir,'Download'))
                if dl_fail_count==5:
                    dl_complete = True
                print('Download sequence failed; reattemping')

        # Check which bands have data
        bands_dict_gal = {}
        for band in bands_dict.keys():
            for raw_file in os.listdir(os.path.join(gal_dir,'Download')):
                if bands_dict[band]['band_short']+'-int.fits' in raw_file:
                    bands_dict_gal[band] = bands_dict[band]

                # Uncompress data (otherwise Montage will fill /tmp/ directory by uncompressing things itself)
                if '.fits.gz' in raw_file:
                    os.system('gzip -d '+os.path.join(gal_dir,'Download',raw_file))

        # Record null file for any band that doens't have data
        for band in bands_dict.keys():
            if band not in bands_dict_gal.keys():
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_GALEX_'+band+'.null'))



        # Loop over bands, to conduct mosaicing
        for band in bands_dict_gal.keys():
            band_dict = bands_dict_gal[band]

             # Declare directories
            id_string = name+'_GALEX_'+band_dict['band_long']

            # Create handling directories for stages of processing
            sub_dir_list = [band+'_Raw', band+'_Diffs_Temp', band+'_Backsub_Temp', band+'_Reproject_Temp', band+'_Convolve_Temp']
            for sub_dir in sub_dir_list:
                if os.path.exists(os.path.join(gal_dir, sub_dir)):
                    shutil.rmtree(os.path.join(gal_dir, sub_dir))
                os.mkdir(os.path.join(gal_dir, sub_dir))

            # Create handling directory for SWarp separately, to stop .NFS locks being an issue
            if not os.path.exists(os.path.join(os.path.join(temp_dir,'SWarp_Temp'))):
                os.mkdir(os.path.join(temp_dir,'SWarp_Temp'))
            swarp_dir = os.path.join(temp_dir,'SWarp_Temp',name+'_'+band+'_'+str(time.time()).replace('.','-'))
            os.mkdir(swarp_dir)

            # Copy maps from relevent band from download directory to raw directory
            for dl_file in os.listdir(os.path.join(gal_dir, 'Download')):
                if '-'+band_dict['band_short']+'-' in dl_file:
                    shutil.copy2(os.path.join(gal_dir, 'Download', dl_file), os.path.join(gal_dir, band+'_Raw'))

            # Ensure that at least one of the raw GALEX int tiles has actual flux coverage at location of source
            coverage = False
            raw_files = []
            [raw_files.append(raw_file) for raw_file in os.listdir(os.path.join(gal_dir, band+'_Raw')) if '-int.fits' in raw_file]
            for raw_file in raw_files:

                # Read in map
                in_image, in_header = astropy.io.fits.getdata(os.path.join(gal_dir, band+'_Raw', raw_file), header=True)

                # Locate pixel coords
                in_wcs = astropy.wcs.WCS(in_header)
                location_pix = in_wcs.wcs_world2pix( np.array([[ np.float(ra), np.float(dec) ]]), 0 )[0]
                pix_i, pix_j = np.int(np.round(location_pix[1])), np.int(np.round(location_pix[0]))

                # Evalulate coverage at location
                if True in [ coord<=0 for coord in [ pix_i-10, pix_i+11, pix_j-10, pix_j+11 ] ]:
                    continue
                try:
                    image_slice = in_image[pix_i-10:pix_i+11, pix_j-10:pix_j+11]
                except:
                    continue
                if np.where(image_slice>0)[0].shape[0]>0:
                    coverage = True

            # If no coverage, just record null file
            if not coverage:
                print('No GALEX '+band_dict['band_long']+' coverage for '+name)
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_GALEX_'+bands_dict[band]['band_long']+'.null'))
                gc.collect()
            elif coverage:

                # Loop over raw int tiles, creating exposure maps, and cleaning images to remove null pixels (also, creating convolved maps for later background fitting)
                print('Cleaning '+str(len(raw_files))+' raw maps for '+id_string)
                pool = mp.Pool(processes=int(0.75*mp.cpu_count()))
                for raw_file in raw_files:
                    pool.apply_async(GALEX_Clean, args=(raw_file, os.path.join(gal_dir,band+'_Raw'), os.path.join(gal_dir,band+'_Reproject_Temp'), os.path.join(gal_dir,band+'_Convolve_Temp'), band_dict,))
                    #GALEX_Clean(raw_file, os.path.join(gal_dir,band+'_Raw'), os.path.join(gal_dir,band+'_Reproject_Temp'), os.path.join(gal_dir,band+'_Convolve_Temp'), band_dict)
                pool.close()
                pool.join()

                # Create FITS header
                pix_width_arcsec = 2.5
                out_hdr = ChrisFuncs.Fits.FitsHeader(ra, dec, width, pix_width_arcsec)

                # Reproject images (dooing individually, instead of using mProjExec, to avoid Montage arbitrarily skipping some images)
                print('Reprojecting '+id_string+' maps')
                os.chdir(os.path.join(gal_dir, band+'_Reproject_Temp'))
                joblib.Parallel( n_jobs=int(0.75*mp.cpu_count()) )\
                               ( joblib.delayed( GALEX_Reproject )\
                               ( id_string, band, out_hdr, list_file, gal_dir, swarp_dir )\
                               for list_file in os.listdir(os.path.join(gal_dir, band+'_Reproject_Temp')) )
                """for list_file in os.listdir(os.path.join(gal_dir, band+'_Reproject_Temp')):
                    GALEX_Reproject( id_string, band, out_hdr, list_file, gal_dir, swarp_dir )"""

                # Rename reprojected files for SWarp
                for list_file in os.listdir(swarp_dir):
                    if '_area.fits' in list_file:
                        os.remove(os.path.join(swarp_dir,list_file))
                    elif 'hdu0_' in list_file:
                        os.rename(os.path.join(swarp_dir,list_file), os.path.join(swarp_dir,list_file.replace('hdu0_','')))
                for list_file in os.listdir(swarp_dir):
                    if '.fits.gz.fits' in list_file:
                        raise Exception('Compressed files found in SWarp directory; something has gone wrong')

                # If more than one image file, commence background-matching
                mosaic_count = len(raw_files)
                if mosaic_count>1:
                    print('Matching background of '+id_string+' maps')
                    ChrisFuncs.Coadd.LevelFITS(swarp_dir, '-int.fits', convfile_dir=os.path.join(gal_dir,band+'_Convolve_Temp'))

                # Use SWarp to co-add images weighted by their error maps
                print('Mosaicing '+id_string+' maps')
                image_width_pixels = str(int((float(width)*3600.)/pix_width_arcsec))
                if swarp_bg == False:
                    swarp_bg = 'N'
                elif swarp_bg == True:
                    swarp_bg = 'Y'
                os.chdir(swarp_dir)
                """subprocess.Popen(['swarp','-IMAGEOUT_NAME', id_string+'_SWarp.fits',
                                  '-WEIGHT_SUFFIX', '.wgt.fits',
                                  '-CENTER_TYPE', 'MANUAL',
                                  '-CENTER', str(ra)+','+str(dec),
                                  '-COMBINE_TYPE', 'WEIGHTED',
                                  '-COMBINE_BUFSIZE', '2048',
                                  '-IMAGE_SIZE ', image_width_pixels+','+image_width_pixels,
                                  '-MEM_MAX', '4096',
                                  '-NTHREADS', str(int(0.75*mp.cpu_count())),
                                  '-RESCALE_WEIGHTS', 'N',
                                  '-RESAMPLE', 'N',
                                  '-SUBTRACT_BACK', swarp_bg,
                                  '-VERBOSE_TYPE', 'QUIET',
                                  '-VMEM_MAX', '4095',
                                  '-WEIGHT_TYPE', 'MAP_WEIGHT'])"""
                os.system('swarp *int.fits -IMAGEOUT_NAME '+id_string+'_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER '+str(ra)+','+str(dec)+' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE '+image_width_pixels+','+image_width_pixels+' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N -RESAMPLE N -SUBTRACT_BACK '+swarp_bg+' -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT')

                # Remove null values, correct for pixel rescaling, and save finalised map to output directory
                in_image, in_header = astropy.io.fits.getdata(os.path.join(swarp_dir, id_string+'_SWarp.fits'), header=True)
                out_image = in_image.copy()
                out_image[np.where(out_image==0)] = np.NaN
                out_image[np.where(out_image<-1E3)] = np.NaN
                out_image[np.where(out_image<=1E-8)] = 0
                out_image *= pix_width_arcsec**2.0 / 1.5**2.0
                astropy.io.fits.writeto(os.path.join(temp_dir, id_string+'.fits'), out_image, header=in_header, overwrite=True)

                # Clean up
                print('Completed mosaicing of '+id_string)
                shutil.rmtree(swarp_dir, ignore_errors=True)
                gc.collect()

        # In parallel, generate final standardised maps for each band
        pool = mp.Pool(processes=2)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( GALEX_Generator, args=(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails,) )
            GALEX_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails)
        pool.close()
        pool.join()

        # Clean memory, and return timings (if more than one target being processed)
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        if len(name) > 1:
            print('Estimated time until GALEX data completed for all targets: '+time_est)

        # Tidy up (best as we can)
        gc.collect()
        try:
            shutil.rmtree(temp_dir)
        except:
            ChrisFuncs.RemoveCrawl(temp_dir)
            print('Unable to fully tidy up temporary directory; probably due to NFS locks on network drive')

    # Report completion
    print('Total time elapsed: '+str((time.time()-time_list[0])/3600.0)+' hours')
    print('All available GALEX imagery acquired for all targets')





# Define function to clean GALEX tiles and create exposure maps
def GALEX_Clean(raw_file, raw_dir, reproj_dir, conv_dir, band_dict):
    print('Cleaning map '+raw_file)

    # Read in image
    in_image, in_header = astropy.io.fits.getdata(os.path.join(raw_dir, raw_file), header=True)
    out_image = in_image.copy()

    # Load and align response map
    rr_image = astropy.io.fits.getdata(os.path.join(raw_dir, raw_file.replace('-int.fits','-rr.fits')))
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = scipy.ndimage.interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean int image using response map
    out_image[ np.where( rr_image<=1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_image = astropy.io.fits.getdata(os.path.join(raw_dir, raw_file.replace('-int.fits','-skybg.fits')))
    bg_zoom = np.float(out_image.shape[0]) / np.float(bg_image.shape[0])
    bg_image = scipy.ndimage.interpolation.zoom(bg_image, bg_zoom, order=0)

    # Clean int image using sky background map
    out_image[ np.where( bg_image<=1E-10 ) ] = np.NaN

    # Set all remaining, and hence "true", zero pixels to be ever-so-slighly non-zero (for ease of later processing)
    out_image += 1E-8

    # Find centre of coverage area
    cov_i = ((np.where( np.isnan(out_image)==False ))[0])
    cov_j = ((np.where( np.isnan(out_image)==False ))[1])
    cov_ellipse = ChrisFuncs.EllipseFit(cov_i, cov_j)
    cov_centre = cov_ellipse[0]
    cov_centre_i, cov_centre_j = cov_centre[0], cov_centre[1]

    # Set all pixels more than 35 arcmin (1400 pizels) from centre to be NaN, as these are typically low-quality
    cov_trim_mask = ChrisFuncs.EllipseMask(out_image, 1400, 1.0, 0.0, cov_centre_i, cov_centre_j)
    out_image[ np.where(cov_trim_mask==0) ] = np.NaN

    # Save cleaned image
    astropy.io.fits.writeto(os.path.join(reproj_dir, raw_file), out_image, header=in_header, overwrite=True )

    # Create convolved version of map, for later use in background-matching
    kernel_width = 100
    """conv_image = scipy.ndimage.filters.gaussian_filter(in_image, kernel_width)
    conv_image[ np.where( np.isnan(out_image)==True ) ] = np.NaN"""
    kernel = astropy.convolution.kernels.Tophat2DKernel(kernel_width)
    conv_image = astropy.convolution.convolve_fft(out_image, kernel, nan_treatment='fill', fill_value=0.0, normalize_kernel=True, allow_huge=True)#, interpolate_nan=True, normalize_kernel=True)
    astropy.io.fits.writeto(os.path.join(conv_dir, raw_file), conv_image, header=in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[np.where(np.isnan(out_image) == False)] = (float(in_header['exptime']))**0.5
    astropy.io.fits.writeto(os.path.join(reproj_dir, raw_file.replace('.fits','.wgt.fits')), exp_image, header=in_header)





# Define function to reproject image and weight files for given source a new projection defined in a header file
def GALEX_Reproject(id_string, band, out_hdr, list_file, gal_dir, swarp_dir):

    # To make sure we go in order and catch problems, we skip weight files in the first instance
    if '.wgt.fits' in list_file:
        return

    # First, reproject the actual integration map; if that fails, just move on with life
    print('Reprojecting map '+list_file)
    os.chdir(os.path.join(gal_dir, band+'_Reproject_Temp'))
    try:
        out_img = reproject.reproject_exact(os.path.join(gal_dir,band+'_Reproject_Temp',list_file), out_hdr, parallel=False)[0]
        astropy.io.fits.writeto(os.path.join(swarp_dir,list_file), data=out_img, header=out_hdr)
    except:
        return

    # Next, reproject the weight map; if that fails, delete the corresponding integration file
    try:
        out_img = reproject.reproject_exact(os.path.join(gal_dir,band+'_Reproject_Temp',list_file.replace('.fits','.wgt.fits')), out_hdr, parallel=False)[0]
        astropy.io.fits.writeto(os.path.join(swarp_dir,list_file.replace('.fits','.wgt.fits')), data=out_img, header=out_hdr)
    except:
        if ('-int.fits' in list_file) and (os.path.exists(os.path.join(gal_dir,band+'_Reproject_Temp',list_file.replace('.fits','.wgt.fits')))):
            os.remove(os.path.join(gal_dir,band+'_Reproject_Temp',list_file.replace('.fits','.wgt.fits')))
        if ('.wgt.fits' in list_file) and (os.path.exists(os.path.join(gal_dir,band+'_Reproject_Temp',list_file.replace('.wgt.fits','.fits')))):
            os.remove(os.path.join(gal_dir,band+'_Reproject_Temp',list_file.replace('.wgt.fits','.fits')))



# Define function to finalise GALEX image of a given source in a given band
def GALEX_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails):
    band = band_dict['band_long']
    wavelength = band_dict['wavelength']
    print('Generating final standardised map of GALEX '+band+' data for '+name)

    # If null file exists for this target in this band, copy it to final output directory
    if os.path.exists(os.path.join(temp_dir,'.'+name+'_GALEX_'+band+'.null')):
        shutil.copy(os.path.join(temp_dir,'.'+name+'_GALEX_'+band+'.null'),
                  os.path.join(out_dir,'.'+name+'_GALEX_'+band+'.null'))
    else:

        # Read in map
        in_img, in_hdr = astropy.io.fits.getdata(os.path.join(temp_dir,name+'_GALEX_'+band+'.fits'), header=True)
        in_wcs = astropy.wcs.WCS(in_hdr)
        in_pix_width_arcsec = 3600.0 * astropy.wcs.utils.proj_plane_pixel_scales(in_wcs).mean()
        out_img = in_img.copy()

        # Calculate conversion factor to get from GALEX native units (counts pix^-1 sec^-1) to Janskies (doing it this way stops inf values appearing)
        unit_factor = ChrisFuncs.ABMagsToJy(band_dict['zero_point'] - ( 2.5 * np.log10(1.0) ))

        # Use conversion factor to convert each pixel value to janskys
        out_img *= unit_factor

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
        out_hdr.set('TELESCOP', 'GALEX', 'Telescope that made this observation')
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
        cutout_wcs.wcs.cdelt = [np.mean([in_hdr['CD1_1'],in_hdr['CD1_2']]), np.mean([in_hdr['CD2_1'],in_hdr['CD2_2']])]
        cutout_wcs.wcs.crval = [in_hdr['CRVAL1'], in_hdr['CRVAL2']]
        cutout_wcs.wcs.ctype = [in_hdr['CTYPE1'], in_hdr['CTYPE2']]
        cutout_wcs_header = cutout_wcs.to_header()
        for card in cutout_wcs_header.cards:
            out_hdr.append(card)

        # Write output FITS file
        astropy.io.fits.writeto(os.path.join(out_dir,name+'_GALEX_'+band+'.fits'), data=out_img, header=out_hdr, overwrite=True)

        # Make thumbnail image of cutout
        if thumbnails:
            try:
                thumb_out = aplpy.FITSFigure(out_dir+name+'_GALEX_'+band+'.fits')
                thumb_out.show_colorscale(cmap='gist_heat', stretch='arcsinh')
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=500, lw=2.5, edgecolor='#01DF3A')
                thumb_out.save(os.path.join(out_dir,name+'_GALEX_'+band+'.jpg'), dpi=125)
                thumb_out.close()
            except:
                print('Failed making thumbnail for '+name)
                pdb.set_trace()

        # Clean memory before finishing
        gc.collect()



# Define function to wget GALEX tiles
def GALEX_wget(tile_url, tile_filename):
    if os.path.exists(tile_filename):
        os.remove(tile_filename)
    success = False
    while success==False:
        try:
            wget.download(tile_url, out=tile_filename)
            success = True
        except:
            print('Failure! Retrying acquistion of '+tile_url)
            time.sleep(0.1)
            success = False



# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")

