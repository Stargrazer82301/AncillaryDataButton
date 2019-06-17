# Identify location
import socket
location = socket.gethostname()
if location == 'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'Hobbitslayer':
    dropbox = 'C:\\Users\\spx7cjc\\Dropbox\\'
if location in ['saruman','scoobydoo','caroline','herschel01','gandalf','ghost']:
    dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Import smorgasbord
import os
import sys
import multiprocessing as mp
import numpy as np
import scipy.ndimage
import astropy.io.fits
import montage_wrapper.commands
import lmfit
import shutil
import signal
import gc
import wget
import pdb
import time
import requests
import xml
import ChrisFuncs





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget GALEX tiles
def GALEX_wget(tile_url, tile_filename, verbose=False):
    if verbose: print 'Acquiring '+tile_url
    if os.path.exists(tile_filename):
        os.remove(tile_filename)
    success = False
    while success==False:
        try:
            wget.download(tile_url, out=tile_filename)
            if verbose: print 'Successful acquisition of '+tile_url
            success = True
        except:
            print 'Failure! Retrying acquistion of '+tile_url
            time.sleep(0.1)
            success = False



# Define function to set a set of maps to the same level
def GALEX_Zero(fitsfile_dir, convfile_dir, target_suffix):

    # Define sill ysubfunction that fits flat plane to image, to find level
    def GALEX_Level_Chisq(level_params, image):
        level = level_params['level'].value
        chi = image - level
        chisq = chi**2.0
        return chisq

    # Make list of files in target directory that have target suffix
    allfile_list = os.listdir(fitsfile_dir)
    fitsfile_list = []
    for allfile in allfile_list:
        if target_suffix in allfile:
            fitsfile_list.append(allfile)

    # Loop over each file
    for i in range(0, len(fitsfile_list)):
        print 'Matching backgorund of map '+fitsfile_list[i]

        # Read in corresponding map from directory containing convolved images
        image_conv = astropy.io.fits.getdata(os.path.join(fitsfile_dir,fitsfile_list[i]))

        # Fit to level of image; save if first image, otherwise calculate appropriate offset
        level_params = lmfit.Parameters()
        level_params.add('level', value=np.nanmedian(image_conv), vary=True)
        image_conv_clipped = ChrisFuncs.SigmaClip(image_conv, tolerance=0.005, median=False, sigma_thresh=3.0)[2]
        level_result = lmfit.minimize(GALEX_Level_Chisq, level_params, args=(image_conv_clipped.flatten(),))
        level = level_result.params['level'].value
        if i==0:
            level_ref = level
            continue
        average_offset = level_ref - level

        """
        # Save floor and peak values
        floor_value = np.nanmin(image_conv)
        peak_value = ChrisFuncs.SigmaClip( image_conv, tolerance=0.00025, median=False, sigma_thresh=3.0)[1]
        floor_value_list.append(floor_value)
        peak_value_list.append(peak_value)
        if i==0:
            floor_value_ref = floor_value
            peak_value_ref = peak_value
            continue

        # Calculate offsets
        floor_offset = floor_value_ref - floor_value
        peak_offset = peak_value_ref - peak_value
        average_offset = peak_offset#np.mean([ floor_offset, peak_offset ])
        """
        # Read in unconvolved file, and apply offset
        image_in, header_in = astropy.io.fits.getdata(os.path.join(fitsfile_dir,fitsfile_list[i]), header=True)
        image_out = image_in + average_offset

        # Save corrected file
        save_success = False
        save_fail = 0
        while not save_success:
            try:
                astropy.io.fits.writeto(os.path.join(fitsfile_dir,fitsfile_list[i]), image_out, header=header_in, overwrite=True)
                save_success = True
            except:
                print 'Unable to save corrected file '+fitsfile_list[i]+'; reattempting'
                save_fail += 1
                time.sleep(10)
                if save_fail>=5:
                    pdb.set_trace()



# Define function to clean GALEX tiles and create exposure maps
def GALEX_Clean(raw_file, raw_dir, reproj_dir, conv_dir, band_dict):
    print 'Cleaning map '+raw_file

    # Read in image
    in_image, in_header = astropy.io.fits.getdata( os.path.join( raw_dir, raw_file ), header=True )
    out_image = in_image.copy()

    # Load and align response map
    rr_image = astropy.io.fits.getdata( os.path.join( raw_dir, raw_file.replace('-int.fits','-rr.fits') ) )
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = scipy.ndimage.interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean int image using response map
    out_image[ np.where( rr_image<=1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_image = astropy.io.fits.getdata( os.path.join( raw_dir, raw_file.replace('-int.fits','-skybg.fits') ) )
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
    astropy.io.fits.writeto( os.path.join( raw_dir, raw_file ), out_image, header=in_header, overwrite=True )

    # Create convolved version of map, for later use in background-matching
    kernel_width = 20
    """conv_image = scipy.ndimage.filters.gaussian_filter(in_image, kernel_width)
    conv_image[ np.where( np.isnan(out_image)==True ) ] = np.NaN"""
    kernel = astropy.convolution.kernels.Tophat2DKernel(kernel_width)
    conv_image = astropy.convolution.convolve_fft(out_image, kernel, nan_treatment='fill', fill_value=0.0, normalize_kernel=True, allow_huge=True)#, interpolate_nan=True, normalize_kernel=True)
    astropy.io.fits.writeto( os.path.join( conv_dir, raw_file ), conv_image, header=in_header )

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5
    astropy.io.fits.writeto( os.path.join( reproj_dir, raw_file.replace('.fits','.wgt.fits') ), exp_image, header=in_header )



# Define function to mosaic together GALEX tiles of a given source
def GALEX_Montage(name, ra, dec, width, band_dict, gal_dir, out_dir):
    print 'Commencing mosaicing for '+name

    # Declare directories
    id_string = name+'_GALEX_'+band_dict['band_long']

    # Create storage directories for Montage and SWarp (deleting any prior), and set appropriate Python working directory
    sub_dir_list = ['Raw', 'Diffs_Temp', 'Backsub_Temp', 'SWarp_Temp', 'Reproject_Temp', 'Convolve_Temp']
    for sub_dir in sub_dir_list:
        if os.path.exists( os.path.join( gal_dir, sub_dir ) ):
            shutil.rmtree( os.path.join( gal_dir, sub_dir ) )
        os.mkdir( os.path.join( gal_dir, sub_dir ) )

    # Copy maps from relevent band from download directory to raw directory
    for dl_file in os.listdir( os.path.join( gal_dir, 'Download' ) ):
        if '-'+band_dict['band_short']+'-' in dl_file:
            shutil.copy2( os.path.join( gal_dir, 'Download', dl_file ), os.path.join( gal_dir, 'Raw' ) )

    # Ensure that at least one of the raw GALEX int tiles has actual flux coverage at location of source
    coverage = False
    raw_files = []
    [ raw_files.append(raw_file) for raw_file in os.listdir( os.path.join( gal_dir, 'Raw') ) if '-int.fits' in raw_file ]
    for raw_file in raw_files:

        # Read in map
        in_image, in_header = astropy.io.fits.getdata( os.path.join( gal_dir, 'Raw', raw_file), header=True )

        # Locate pixel coords
        in_wcs = astropy.wcs.WCS(in_header)
        location_pix = in_wcs.wcs_world2pix( np.array([[ np.float(ra), np.float(dec) ]]), 0 )[0]
        pix_i, pix_j = np.int(np.round(location_pix[1])), np.int(np.round(location_pix[0]))

        # Evalulate coverage at location, and proceed accordingly
        if True in [ coord<=0 for coord in [ pix_i-10, pix_i+11, pix_j-10, pix_j+11 ] ]:
            continue
        try:
            image_slice = in_image[pix_i-10:pix_i+11, pix_j-10:pix_j+11]
        except:
            continue
        if np.where(image_slice>0)[0].shape[0]>0:
            coverage = True
    if not coverage:
        print 'No GALEX '+band_dict['band_long']+' coverage for '+name
        gc.collect()
    elif coverage:

        # Loop over raw int tiles, creating exposure maps, and cleaning images to remove null pixels (also, creating convolved maps for later background fitting)
        print 'Cleaning '+str(len(raw_files))+' raw maps for '+id_string
        pool = mp.Pool(processes=7)
        for raw_file in raw_files:
            pool.apply_async(GALEX_Clean, args=(raw_file, os.path.join(gal_dir,'Raw'), os.path.join(gal_dir,'Reproject_Temp'), os.path.join(gal_dir,'Convolve_Temp'), band_dict,))
            #GALEX_Clean(raw_file, os.path.join(gal_dir,'Raw'), os.path.join(gal_dir,'Reproject_Temp'), os.path.join(gal_dir,'Convolve_Temp'), band_dict)
        pool.close()
        pool.join()

        # Create Montage FITS header
        location_string = str(ra)+' '+str(dec)
        pix_size = 2.5
        montage_wrapper.commands.mHdr(location_string, width, os.path.join( gal_dir, id_string+'_HDR' ), pix_size=pix_size)

        # Count image files, and move to reprojection directory
        mosaic_count = len(raw_files)
        for raw_file in raw_files:
            shutil.move( os.path.join( gal_dir, 'Raw', raw_file ), os.path.join( gal_dir, 'Reproject_Temp' ) )

        # If more than one image file, commence background-matching
        if mosaic_count>1:
            print 'Matching background of '+id_string+' maps'
            GALEX_Zero( os.path.join( gal_dir, 'Reproject_Temp' ), os.path.join( gal_dir, 'Convolve_Temp' ), '-int.fits' )

        # Reproject image and weight prior to coaddition
        os.chdir( os.path.join( gal_dir, 'Reproject_Temp' ) )
        montage_wrapper.commands.mImgtbl( os.path.join(gal_dir,'Reproject_Temp'),  os.path.join(gal_dir,band+'_Image_Metadata_Table.dat'), corners=True)
        montage_wrapper.commands.mProjExec( os.path.join(gal_dir,band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,id_string+'_HDR'), os.path.join(gal_dir,'SWarp_Temp'), os.path.join(gal_dir,id_string+'_Proj_Stats.txt'), raw_dir=os.path.join(gal_dir,'Reproject_Temp'), debug=False, exact=True, whole=False)

        # Rename reprojected files for SWarp
        for listfile in os.listdir( os.path.join(gal_dir,'SWarp_Temp') ):
            if '_area.fits' in listfile:
                os.remove( os.path.join(gal_dir,'SWarp_Temp',listfile) )
            elif 'hdu0_' in listfile:
                os.rename( os.path.join(gal_dir,'SWarp_Temp',listfile), os.path.join(gal_dir,'SWarp_Temp',listfile.replace('hdu0_','')) )
        for listfile in os.listdir( os.path.join(gal_dir,'SWarp_Temp') ):
            if '.fits.gz.fits' in listfile:
                os.rename( os.path.join(gal_dir,'SWarp_Temp',listfile), os.path.join(gal_dir,'SWarp_Temp',listfile.replace('.fits.gz.fits','.fits')) )

        # Use SWarp to co-add images weighted by their error maps
        print 'Co-adding '+id_string+' maps'
        image_width_pixels = str(int((float(width)*3600.)/pix_size))
        os.chdir( os.path.join( gal_dir, 'SWarp_Temp' ) )
        os.system('swarp *int.fits -IMAGEOUT_NAME '+id_string+'_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER '+str(ra)+','+str(dec)+' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE '+image_width_pixels+','+image_width_pixels+' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT')

        # Remove null values, correct for pixel rescaling, and save finalised map to output directory
        in_image, in_header = astropy.io.fits.getdata( os.path.join( gal_dir, 'SWarp_Temp', id_string+'_SWarp.fits' ), header=True )
        out_image = in_image.copy()
        out_image[ np.where( out_image==0 ) ] = np.NaN
        out_image[ np.where( out_image<-1E3 ) ] = np.NaN
        out_image[ np.where( out_image<=1E-8 ) ] = 0
        out_image *= pix_size**2.0 / 1.5**2.0
        astropy.io.fits.writeto( os.path.join( out_dir, id_string+'.fits' ), out_image, header=in_header, overwrite=True )

        # Clean up
        print 'Completed Montaging and SWarping of '+id_string
        gc.collect()





# Commence main task
if __name__ == "__main__":

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/GALEX/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/GALEX/Mosaics/'

    # Read in source catalogue
    NESS_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = NESS_cat['name']

    # State band information
    bands_dict = {'FUV':{'band_short':'fd','band_long':'FUV'},
                  'NUV':{'band_short':'nd','band_long':'NUV'}}

    # Provide information about different map types
    map_types = ['Raw','Background','Exposure','Response']#,'Flags']
    map_suffixes = ['-int.fits.gz','-skybg.fits.gz','-exp.fits.gz','-rr.fits.gz']#,'-flags.fits.gz']

    # Delete any pre-existing temporary files
    if os.path.exists(in_dir):
        shutil.rmtree(in_dir)
    os.mkdir(in_dir)

    # Read in list of already-processed sources, and identify sources not yet processed
    already_file = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/GALEX/GALEX_Already_Processed_List.dat'
    if not os.path.exists(already_file):
        open(already_file,'a')
    alrady_processed = np.genfromtxt(already_file, dtype=('S50')).tolist()
    remaining_list = []
    for i in range(0, name_list.shape[0]):
        already_done = 0
        name = name_list[i]
        if name not in alrady_processed:
            remaining_list.append(i)
    NESS_cat = NESS_cat[remaining_list]
    name_list = NESS_cat['name']
    ra_list = NESS_cat['ra']
    dec_list = NESS_cat['dec']

    # Register signal function handler, for dealing with timeouts
    signal.signal(signal.SIGALRM, Handler)

    # Record time taken
    time_list = [time.time()]

    # Loop over each source
    for i in range(0, NESS_cat.shape[0]):
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        width = 0.25
        print 'Processing target '+name

        # Try to catch chaos in the act, after multiple retries
        target_complete = False
        target_fail_count = 0
        while not target_complete:
            try:
                signal.alarm(1800)

                # Create tile processing dirctories
                os.mkdir( os.path.join( in_dir, name ) )
                os.mkdir( os.path.join( in_dir, name, 'Download' ) )

                # Perform query (removing pre-existing query file, if present)
                print 'Querying GALEX server using SIA protocol'
                query_success = False
                while not query_success:
                    try:
                        query_url = 'http://mast.stsci.edu/portal_vo/Mashup/VoQuery.asmx/SiaV1?MISSION=GALEX&POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS&FORMAT=image/fits'
                        query_filename = os.path.join(in_dir,name,str(name)+'.vot')
                        if os.path.exists(query_filename):
                            os.remove(query_filename)
                        wget.download(query_url, out=query_filename)
                        query_success = True
                    except:
                        query_success = False

                # Read query result VOTable
                query_table = astropy.table.Table.read(query_filename)

                # If query finds no matches, continue to next target
                if len(query_table)==0:
                    print 'No GALEX data for '+name
                    alrady_processed_file = open(already_file, 'a')
                    alrady_processed_file.write(name+'\n')
                    alrady_processed_file.close()
                    shutil.rmtree( os.path.join( in_dir, name ) )
                    time_list.append(time.time())
                    target_complete = True
                    continue

                # Grab fits files URLs
                print 'Processing query result table'
                query_urls = []
                for i in range(len(query_table)):
                    query_urls.append(query_table[i]['accessURL'])

                # In parallel, download files to holding directory
                print 'Retrieving GALEX data for '+name
                os.chdir(in_dir)
                dl_pool = mp.Pool(processes=20)
                for j in range(0, len(query_urls)):
                    primary_url = query_urls[j]
                    for map_suffix in map_suffixes:
                        tile_url = primary_url.replace('-int.fits.gz', map_suffix)
                        tile_filename = os.path.join( in_dir, name, 'Download', tile_url.split('/')[-1] )
                        #GALEX_wget( tile_url, os.path.join(in_dir,name,'Download',tile_filename) )
                        dl_pool.apply_async( GALEX_wget, args=( tile_url, os.path.join(in_dir,name,'Download',tile_filename), ) )
                dl_pool.close()
                dl_pool.join()

                # Check which bands have data
                bands_dict_source = {}
                for band in bands_dict.keys():
                    for raw_file in os.listdir( os.path.join( in_dir, name, 'Download' ) ):
                        if bands_dict[band]['band_short']+'-int.fits.gz' in raw_file:
                            bands_dict_source[band] = bands_dict[band]

                # Loop over bands, conducting mosaicing function
                for band in bands_dict_source.keys():
                    GALEX_Montage(name, ra, dec, width, bands_dict_source[band], os.path.join( in_dir, name ), out_dir )

            except:
                print('Failure!')

            # Record that processing of souce has been compelted
            alrady_processed_file = open(already_file, 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()

            # Clean memory, and return timings
            shutil.rmtree( os.path.join( in_dir, name ) )
            gc.collect()
            time_list.append(time.time())
            time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
            time_file = open( os.path.join('/'.join(in_dir.split('/')[:-2]),'Estimated_Completion_Time.txt'), 'w')
            time_file.write(time_est)
            time_file.close()
            target_complete = True
            print 'Estimated completion time: '+time_est


# Jubilate
print 'All done!'



