# Identify location
import socket
location = socket.gethostname()
if location  ==  'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location  ==  'saruman':
    dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Import smorgasbord
import pdb
import os
import sys
import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import astropy.io.votable
import montage_wrapper.commands
import shutil
import signal
import gc
#warnings.filterwarnings("ignore")
import subprocess
import wget
import time
import ChrisFuncs
import ChrisFuncs.Coadd

# Add SWarp directory to path
os.environ['PATH'] = os.environ['PATH'] + ':/home/user/spx7cjc/swarp/bin'





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget and extact SMSS files
def SMSS_wget(data_url, data_filename):
    print 'Acquiring '+data_url
    if os.path.exists(data_filename):
        os.remove(data_filename)
    success = False
    while success == False:
        try:
            wget.download(data_url, out=data_filename)
            print 'Successful acquisition of '+data_url
            success = True
        except:
            print 'Failure! Retrying acquistion of '+data_url
            time.sleep(5.0)
            success = False



# Define function to Montage together contents of folder
def SMSS_Montage(name, ra, dec, pix_width, map_width, band, input_dir, out_dir):
    print 'Montaging '+name+'_SMSS_'+band

    # Make temporary directories
    if os.path.exists(os.path.join(input_dir,'Montage_Temp')): shutil.rmtree(os.path.join(input_dir,'Montage_Temp'))
    os.mkdir(os.path.join(input_dir,'Montage_Temp'))
    if os.path.exists(os.path.join(input_dir,'Montage_Work')): shutil.rmtree(os.path.join(input_dir,'Montage_Work'))
    os.mkdir(os.path.join(input_dir,'Montage_Work'))
    if os.path.exists(os.path.join(input_dir,'Proj_Temp')): shutil.rmtree(os.path.join(input_dir,'Proj_Temp'))
    os.mkdir(os.path.join(input_dir,'Proj_Temp'))
    if os.path.exists(os.path.join(input_dir,'Diff_Temp')): shutil.rmtree(os.path.join(input_dir,'Diff_Temp'))
    os.mkdir(os.path.join(input_dir,'Diff_Temp'))
    if os.path.exists(os.path.join(input_dir,'Corr_Temp')): shutil.rmtree(os.path.join(input_dir,'Corr_Temp'))
    os.mkdir(os.path.join(input_dir,'Corr_Temp'))
    if os.path.exists(os.path.join(input_dir,'SWarp_Temp')): shutil.rmtree(os.path.join(input_dir,'SWarp_Temp'))
    os.mkdir(os.path.join(input_dir,'SWarp_Temp'))
    location_string = str(ra)+' '+str(dec)
    os.chdir(input_dir)

    # Regrid images to same projection
    montage_wrapper.commands.mHdr(location_string, map_width, os.path.join(input_dir,str(name)+'_HDR'), pix_size=0.4)
    montage_wrapper.commands.mImgtbl(os.path.join(input_dir), os.path.join(input_dir,band+'_Image_Metadata_Table.dat'), corners=True)
    montage_wrapper.commands.mProjExec(os.path.join(input_dir,band+'_Image_Metadata_Table.dat'), os.path.join(input_dir,str(name)+'_HDR'), os.path.join(input_dir,'Proj_Temp'), os.path.join(input_dir,band+'_Proj_Stats.txt'), raw_dir=os.path.join(in_dir), debug=False, exact=True, whole=True)
    os.chdir(os.path.join(input_dir,'Proj_Temp'))
    os.system("rename.ul hdu0_IRAS IRAS *IRAS*")
    os.chdir(input_dir)

    # Model and match backgrounds for all images
    montage_wrapper.commands.mOverlaps(os.path.join(input_dir,band+'_Image_Metadata_Table.dat'), os.path.join(input_dir,band+'_Image_Diffs_Table.dat'))
    montage_wrapper.commands.mDiffExec(os.path.join(input_dir,band+'_Image_Diffs_Table.dat'), os.path.join(input_dir,str(name)+'_HDR'), os.path.join(input_dir,'Diff_Temp'), proj_dir=os.path.join(input_dir,'Proj_Temp'), no_area=True)
    montage_wrapper.commands.mFitExec(os.path.join(input_dir,band+'_Image_Diffs_Table.dat'), os.path.join(input_dir,band+'_Image_Fitting_Table.dat'), os.path.join(input_dir,'Diff_Temp'))
    montage_wrapper.commands.mBgModel(os.path.join(input_dir,band+'_Image_Metadata_Table.dat'), os.path.join(input_dir,band+'_Image_Fitting_Table.dat'), os.path.join(input_dir,band+'_Image_Corrections_Table.dat'), level_only=True, n_iter=16384)
    montage_wrapper.commands.mBgExec(os.path.join(input_dir,band+'_Image_Metadata_Table.dat'), os.path.join(input_dir,band+'_Image_Corrections_Table.dat'), os.path.join(input_dir,'Corr_Temp'), proj_dir=os.path.join(input_dir,'Proj_Temp'), no_area=True)

    # Rename corrected files for SWarp to handle
    os.chdir(os.path.join(input_dir,'Corr_Temp'))
    os.system('rename.ul .fits _smss.fits *.fits')

    for in_file in os.listdir(os.path.join(input_dir,'Corr_Temp')):
        if '.fits' not in in_file:
            continue
        ChrisFuncs.FitsCutout(os.path.join(input_dir,'Corr_Temp',in_file), ra, dec, 0.5*map_width*3600.0, pix_width_arcsec=0.5, reproj=True, outfile=os.path.join(input_dir,'SWarp_Temp',in_file), parallel=False, fast=False)

    # Create weight maps using exposure times in headers
    for in_file in os.listdir(os.path.join(input_dir,'SWarp_Temp')):
        if '.fits' not in in_file:
            continue
        in_map, in_header = astropy.io.fits.getdata(os.path.join(input_dir,'SWarp_Temp',in_file), header=True)
        in_exp = astropy.io.fits.getheader(os.path.join(input_dir,in_file.replace('_smss.fits','.fits')))['EXPTIME']
        in_wgt = in_exp**0.5
        out_map = np.zeros(in_map.shape)
        out_map[np.where(np.isnan(in_map))] = np.nan
        out_map[np.where(np.isnan(in_map)==False)] = in_wgt
        astropy.io.fits.writeto(os.path.join(input_dir,'SWarp_Temp',in_file.replace('_smss.fits','_smss.wgt.fits')), data=out_map, header=in_header)

    # Coadd images using weights
    os.chdir(os.path.join(input_dir,'SWarp_Temp'))
    os.system('swarp *_smss.fits -IMAGEOUT_NAME '+name+'_SMSS_'+band+'_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -GAIN_KEYWORD DIESPIZERDIE -RESCALE_WEIGHTS N -SUBTRACT_BACK N -RESAMPLE N -VMEM_MAX 4095 -MEM_MAX 4096 -WEIGHT_TYPE MAP_WEIGHT -NTHREADS 4 -VERBOSE_TYPE QUIET')

    # handle NaN pixels, and write to output directoryu
    swarp_map, swarp_header = astropy.io.fits.getdata(os.path.join(input_dir,'SWarp_Temp',name+'_SMSS_'+band+'_SWarp.fits'), header=True)
    swarp_map[np.where(swarp_map<-1E25)] = np.nan
    astropy.io.fits.writeto(os.path.join(out_dir,name+'_SMSS_'+band+'.fits'), data=swarp_map, header=swarp_header, overwrite=True)





# Commence main task
if __name__  ==  "__main__":

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SMSS/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SMSS/Mosaics/'

    # Read in source catalogue
    nesscat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = nesscat['name']

    # Read in list of already-Montaged sources
    already_file = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SMSS/SMSS_Already_Processed_List.dat'
    if not os.path.exists(already_file):
        open(already_file,'a')
    already_processed = np.genfromtxt(already_file, dtype=('S50')).tolist()

    # Identify targets not yet processed
    remaining_list = []
    for i in range(0, name_list.shape[0]):
        already_done = 0
        name = name_list[i]
        if name not in already_processed:
            remaining_list.append(i)
    name_list = nesscat['name']
    ra_list = nesscat['ra']
    dec_list = nesscat['dec']

    # State band information
    bands_list = ['u','v','g','r','i','z']
    bands_zp = [25.47, 24.88, 25.68, 25.46, 25.43, 25.47]

    # Register signal function handler, for dealing with timeouts
    signal.signal(signal.SIGALRM, Handler)

    # Record time taken
    time_list = [ time.time() ]

    # Loop over each source
    for i in np.random.permutation(np.array(range(0, nesscat.shape[0]))):
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        time_start = time.time()
        width = 0.25
        print 'Processing source '+name

        # Create data processing dirctory (deleting any prior)
        gal_dir = in_dir+str(name)+'/'
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)

        # If source is in list of already-montaged sources, or too far north for SkyMapper to observe, continue
        if (ra > 15.0) or (name in already_processed):
            if ra > 15.0:
                print name+' too far north to be observed by SkyMapper'
            if name in already_processed:
                print name+' already processed'
            alrady_processed_file = open(already_file, 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()
            shutil.rmtree(gal_dir)
            time_list.append(time.time())
            continue

        # Loop over bands
        bands_data = [False]*len(bands_list)
        for b in range(len(bands_list)):
            band = bands_list[b]
            smss_urls = []
            os.makedirs(os.path.join(gal_dir,band))

            # To get around the SMSS 10' limit, set up a grid of sky positions to query
            if width > 9.0/60.0:
                ra_range = np.linspace(ra-width, ra+width, num=2*np.ceil((width/(9.0/60.0))+1))
                dec_range = np.linspace(dec-width, dec+width, num=2*np.ceil((width/(9.0/60.0))+1))
            else:
                ra_range = np.array([ra])
                dec_range = np.array([dec])

            # Loop over sky positions in grid, to query in term
            ra_dec_grid = np.zeros([len(ra_range),len(dec_range)])
            for i in range(len(ra_range)):
                for j in range(len(dec_range)):
                    ra_temp = ra_range[i]
                    dec_temp = dec_range[j]

                    # Construct URL for this query
                    query_url = 'http://skymappersiap.asvo.nci.org.au/dr1_cutout/query?POS='+str(ra_temp)+','+str(dec_temp)+'&SIZE='+str(10.0/60.0)+'&BAND='+band+'&FORMAT=image/fits&INTERSECT=OVERLAPS'

                    # Perform query, with error handling
                    print 'Querying SMSS SIAP server at position '+str(ra_temp)+' '+str(dec_temp)
                    query_success = False
                    query_fail_count = 0
                    while query_success == False:
                        if query_fail_count>=5:
                            break
                        try:
                            query_filename = os.path.join(in_dir,name,str(name)+'.vot')
                            if os.path.exists(query_filename):
                                os.remove(query_filename)
                            wget.download(query_url, out=query_filename)

                            # Read query result VOTable
                            query_output = astropy.io.votable.parse_single_table(query_filename)
                            query_table = query_output.array

                            # Loop over results, recording data download URLs
                            for query_row in query_table:
                                smss_urls.append(query_row['get_image'])

                            # Record success, or handle query failure
                            query_success = True
                            ra_dec_grid[i,j] = 1
                        except:
                            print 'SIAP query failed; reattempting'
                            query_fail_count += 1
                            time.sleep(5)
                    if not os.path.exists(query_filename):
                        query_fail_count += 1
                        query_success = False

            # If queries failed, or no URLs, continue to next band
            smss_urls = list(set(smss_urls))
            if (ra_dec_grid.sum() == 0) or (len(smss_urls) == 0):
                print 'All SIAP queries failed; progressing to next target'
                continue
            bands_data[b] = True

            # In parallel, download files
            os.chdir(os.path.join(gal_dir,band))
            dl_pool = mp.Pool(processes=20)
            for j in range(0, len(smss_urls)):
                data_url = smss_urls[j]
                data_filename = os.path.join(os.path.join(gal_dir,band),name+'_'+band+'_'+str(j)+'.fits')
                dl_pool.apply_async( SMSS_wget, args=(data_url, data_filename,) )
                #SMSS_wget(data_url, data_filename)
            dl_pool.close()
            dl_pool.join()

        # If query finds no matches, continue to next target
        if max(bands_data) == False:
            print 'No SMSS data for '+name
            alrady_processed_file = open(already_file, 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()
            shutil.rmtree( os.path.join( in_dir, name ) )
            time_list.append(time.time())
            continue



        # Loop over bands
        for b in range(len(bands_list)):
            band = bands_list[b]
            band_listmaps = os.listdir(os.path.join(gal_dir,band))
            if len(band_listmaps) == 0:
                continue

            # Loop over maps, first checking they are valid
            for raw_fits in band_listmaps:
                try:
                    test_image, test_header = astropy.io.fits.getdata(os.path.join(gal_dir,band,raw_fits), header=True)
                    print 'Preprocessing file: '+raw_fits
                except Exception, exception_msg:
                    print exception_msg
                    print 'Invalid file: '+raw_fits
                    continue

                # Load file, and retain original header for later use
                in_image, in_header = astropy.io.fits.getdata(os.path.join(gal_dir,band,raw_fits), header=True)
                raw_wcs = astropy.wcs.WCS(in_header)

                # Use header information and Vega-to-AB conversions to work out true zero-point magnitude, as per prescription here: https://www.eso.org/sci/observing/phase3/data_releases/vvv_dr4.pdf
                header_exptime = in_header['EXPTIME']
                band_zp = bands_zp[b]

                # Set all pixels to have units of AB mag, using the calculated zero-point magnitude
                out_image = band_zp - ( 2.5 * np.log10(in_image) )

                # Convert pixel units to Jy
                out_image = ChrisFuncs.ABMagsToJy(out_image)

                # Save converted image
                astropy.io.fits.writeto(os.path.join(gal_dir,band,raw_fits), out_image, header=in_header, clobber=True)



        # In parallel, Montage together each band's input tiles, whilst dealing with timeouts
        complete = False
        fail_counter = 0
        while not complete:
            if fail_counter >= 3:
                continue
#            try:
#            signal.alarm(7200)
            pool = mp.Pool(processes=6)
            for b in range(0, len(bands_list)):
                if bands_data[b]==False:
                    continue
                band = bands_list[b]
                pool.apply_async( SMSS_Montage, args=(name, ra, dec, 0.5, width, band, os.path.join( in_dir, name, band ), out_dir,) )
                #SMSS_Montage(name, ra, dec, 0.5, width, band, os.path.join( in_dir, name, band ), out_dir)#3600.0*pix_min[b]
            pool.close()
            pool.join()
            shutil.rmtree( os.path.join( in_dir, name ) )
            complete = True
#            except Exception, exception_msg:
#                fail_counter += 1
#                gc.collect()
#                for band in bands_list:
#                    input_dir = os.path.join(in_dir,name)+band+'/'
#                    if os.path.exists(input_dir+'/Montage_Temp'):
#                        shutil.rmtree(input_dir+'/Montage_Temp')
#                print 'Mosaicing failure!'
#                print exception_msg

        # Record that processing of souce has been compelted
        alrady_processed_file = open(already_file, 'a')
        alrady_processed_file.write(name+'\n')
        alrady_processed_file.close()

        # Clean memory, and return timings
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        time_file = open( os.path.join('/'.join(in_dir.split('/')[:-2]),'Estimated_Completion_Time.txt'), 'w')
        time_file.write(time_est)
        time_file.close()
        print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'
