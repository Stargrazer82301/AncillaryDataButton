# Identify location
import socket
location = socket.gethostname()
if location == 'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'Hobbitslayer':
    dropbox = 'C:\\Users\\spx7cjc\\Dropbox\\'
if location == 'saruman':
    dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Import smorgasbord
import os
import sys
import multiprocessing as mp
import numpy as np
import astropy.io.votable
from astroquery import sha
import montage_wrapper.commands
import shutil
import signal
import gc
#warnings.filterwarnings("ignore")
from glob import glob
import subprocess
import wget
import pdb
import time
import ChrisFuncs
import ChrisFuncs.Coadd

# Add SWarp directory to path
os.environ['PATH'] = os.environ['PATH'] + ':/home/user/spx7cjc/swarp/bin'





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget and extact Spitzer files
def Spitzer_wget(tile_url, tile_filename):
    print 'Acquiring '+tile_url
    if os.path.exists(tile_filename):
        os.remove(tile_filename)
    success = False
    while success==False:
        try:
            wget.download(tile_url, out=tile_filename)
            print 'Successful acquisition of '+tile_url
            success = True
        except:
            print 'Failure! Retrying acquistion of '+tile_url
            time.sleep(0.1)
            success = False
    os.system('unzip '+tile_filename)



# Defien function to replace null pixels in SWarp outputs with NaNs
def Spitzer_SWarp_NaN(target):
    in_fitsdata = astropy.io.fits.open(target)
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    out_image = in_image.copy()
    out_image[ np.where( out_image<-1E20 ) ] = np.NaN
    out_hdu = astropy.io.fits.PrimaryHDU(data=out_image, header=in_header)
    out_hdulist = astropy.io.fits.HDUList([out_hdu])
    out_hdulist.writeto(target, clobber=True)





# Commence main task
if __name__ == "__main__":

    # Decide what intrument to work for
    instrument = 'IRAC'

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Spitzer/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Spitzer/Mosaics_'+instrument+'/'

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']

    # Read in list of already-Montaged sources
    already_processed_path = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Spitzer/Spitzer_'+instrument+'_Already_Processed_List.dat'
    if not os.path.exists(already_processed_path):
        open(already_processed_path,'a')
    already_processed = np.genfromtxt(already_processed_path, dtype=('S50')).tolist()

    # Identify targets not yet processed
    remaining_list = []
    for i in range(0, name_list.shape[0]):
        already_done = 0
        name = name_list[i]
        if name not in already_processed:
            remaining_list.append(i)
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # State band information
    if instrument=='IRAC':
        bands_dict = {'3.6um':{'instrument':'IRAC','band_long':'3.6','channel':'ch1','pix_size':0.6},
                      '4.5um':{'instrument':'IRAC','band_long':'4.5','channel':'ch2','pix_size':0.6},
                      '5.8um':{'instrument':'IRAC','band_long':'5.8','channel':'ch3','pix_size':0.6},
                      '8.0um':{'instrument':'IRAC','band_long':'8.0','channel':'ch4','pix_size':0.6}}
    elif instrument=='MIPS':
        bands_dict = {'24um':{'instrument':'MIPS','band_long':'24','channel':'ch1','pix_size':2.45},
                      '70um':{'instrument':'MIPS','band_long':'70','channel':'ch2','pix_size':4.0},
                      '160um':{'instrument':'MIPS','band_long':'160','channel':'ch3','pix_size':8.0}}

    # Register signal function handler, for dealing with timeouts
    signal.signal(signal.SIGALRM, Handler)

    # Record time taken
    time_list = [ time.time() ]

    # Loop over each source
    for i in np.random.permutation(range(0, ness_cat.shape[0])):
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        time_start = time.time()
        width = 0.25#
        print 'Processing source '+name

        # Check if source is in list of already-montaged sources
        if name in already_processed:
            print name+' already processed'
            continue

        # Check which, if any, bands already have data
        print 'Checking existing finalised cutouts for matches to current source'
        bands_dict_req = {}
        for band in bands_dict.keys():
            if name+'_Spitzer_'+bands_dict[band]['band_long']+'.fits.gz' not in os.listdir('/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Spitzer/Cutouts/'):
                bands_dict_req[band] = bands_dict[band]

        # Create tile processing dirctories (deleting any prior), and set appropriate Python (ie, Montage) working directory
        gal_dir = in_dir+str(name)+'/'
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)
        os.makedirs(gal_dir+'Errors')
        os.makedirs(gal_dir+'Raw')
        os.chdir(gal_dir+'Raw')

        # Perform query, with error handling
        print 'Querying Spitzer server'
        query_success = False
        query_fail_count = 0
        while query_success==False:
            if query_fail_count>=10:
                break
            try:
                print 'NOTE: Astroquery currently not working with Spitzer; gettings query results using Spitzer API instead'
                """
                query_obj = Spitzer.query(ra=ra, dec=dec, size=width)
                """
                query_url = 'http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/DataService?RA='+str(ra)+'&DEC='+str(dec)+'&SIZE='+str(width)+'&VERB=3&DATASET=ivo%3A%2F%2Firsa.csv%2Fspitzer.level2'
                query_filename = gal_dir+'Spitzer_Query.csv'
                if os.path.exists(query_filename):
                    os.remove(query_filename)
                wget.download(query_url, out=query_filename.replace('.csv','.txt'))
                os.system('stilts  tcopy ifmt=ipac ofmt=csv '+query_filename.replace('.csv','.txt')+' '+query_filename)
                os.remove(query_filename.replace('.csv','.txt'))
                query_success = True
            except:
                print 'Spitzer query failed; reattempting'
                query_fail_count += 1
        if not os.path.exists(query_filename):
            query_success=False
        if query_success==False:
            print 'No Spitzer data for '+name
            already_processed_file = open(already_processed_path, 'a')
            already_processed_file.write(name+'\n')
            already_processed_file.close()
            time_list.append( time.time() )
            shutil.rmtree(gal_dir)
            gc.collect()
            continue

        # Save query result (removing pre-existing query file, if present)
        """
        query_filename = gal_dir+'Spitzer_Query.csv'
        query_obj.write(query_filename, format='csv')
        """

        # Establish if any data was found; if not, skip
        query_in = np.genfromtxt(query_filename, delimiter=',', names=True, dtype=None)
        if query_in.size==0:
            print 'No Spitzer data for '+name
            already_processed_file = open(already_processed_path, 'a')
            already_processed_file.write(name+'\n')
            already_processed_file.close()
            time_list.append( time.time() )
            shutil.rmtree(gal_dir)
            gc.collect()
            continue

        # Record which urls correspond to data in the desired bands (dealing with awkwardness for if there is only 1 entry, or silly massive files)
        Spitzer_urls = []
        Spitzer_bands = []
        if query_in.size==1:
            if query_in['accessWithAnc1Url']!='NONE' and query_in['filesize']<1E9:
                for band in bands_dict_req.keys():
                    if query_in['wavelength']==bands_dict_req[band]['instrument']+' '+band:
                        Spitzer_urls.append(query_in['accessWithAnc1Url'])
                        Spitzer_bands.append(band)
        else:
            for j in range(0, query_in.size):
                if query_in[j]['accessWithAnc1Url']!='NONE' and query_in[j]['filesize']<1E9:
                    for band in bands_dict_req.keys():
                        if query_in[j]['wavelength']==bands_dict_req[band]['instrument']+' '+band:
                            Spitzer_urls.append(query_in[j]['accessWithAnc1Url'])
                            Spitzer_bands.append(band)

        # In parallel, download and extract files
        os.chdir(gal_dir+'Raw')
        dl_pool = mp.Pool(processes=20)
        for j in range(0, len(Spitzer_urls)):
            tile_url = Spitzer_urls[j]
            tile_filename = gal_dir+'Raw/'+name+'_'+Spitzer_bands[j]+'_'+str(j)+'.zip'
            dl_pool.apply_async( Spitzer_wget, args=(tile_url, tile_filename,) )#Spitzer_wget(tile_url, tile_filename)
        dl_pool.close()
        dl_pool.join()
        [ os.remove(dl_zip) for dl_zip in os.listdir(gal_dir+'Raw/') if '.zip' in dl_zip ]

        # Copy files to relevant folders
        for band in bands_dict_req.keys():
            if os.path.exists(gal_dir+band+'/'):
                shutil.rmtree(gal_dir+band+'/')
            if os.path.exists(gal_dir+'Errors/'+band+'/'):
                shutil.rmtree(gal_dir+'Errors/'+band+'/')
            os.makedirs(gal_dir+band+'/')
            os.makedirs(gal_dir+'Errors/'+band+'/')
            channel = bands_dict_req[band]['channel']
            for dl_folder in os.listdir(gal_dir+'Raw/'):
                if os.path.exists(gal_dir+'Raw/'+dl_folder+'/'+channel+'/pbcd'):
                    for dl_file in os.listdir(gal_dir+'Raw/'+dl_folder+'/'+channel+'/pbcd'):
                        if '_maic.fits' in dl_file:
                            shutil.copy2(gal_dir+'Raw/'+dl_folder+'/'+channel+'/pbcd/'+dl_file, gal_dir+band)
                        if '_munc.fits' in dl_file:
                            shutil.copy2(gal_dir+'Raw/'+dl_folder+'/'+channel+'/pbcd/'+dl_file, gal_dir+band)
        shutil.rmtree(gal_dir+'Raw')

        # Check that the retrieved files provide actual coverage of the point in question
        coverage_bands = []
        for band in bands_dict_req.keys():
            montage_wrapper.commands.mImgtbl(gal_dir+band,  gal_dir+band+'/'+band+'_Image_Metadata_Table.dat', corners=True)
            if os.stat(gal_dir+band+'/'+band+'_Image_Metadata_Table.dat').st_size==0:
                continue
            montage_wrapper.commands_extra.mCoverageCheck(gal_dir+band+'/'+band+'_Image_Metadata_Table.dat', gal_dir+band+'/'+band+'_Overlap_Check.dat', mode='point', ra=ra, dec=dec)
            if sum(1 for line in open(gal_dir+band+'/'+band+'_Overlap_Check.dat'))>3:
                coverage_bands.append(band)
        if len(coverage_bands)==0:
            print 'No Spitzer data for '+name
            already_processed_file = open(already_processed_path, 'a')
            already_processed_file.write(name+'\n')
            already_processed_file.close()
            time_list.append( time.time() )
            shutil.rmtree(gal_dir)
            gc.collect()
            continue

        # Loop over each band for coaddition
        for band in coverage_bands:
            print 'Commencing Montaging and SWarping of '+name+'_Spitzer_'+band
            os.chdir(gal_dir+band)
            os.mkdir(gal_dir+band+'/Diffs_Temp')
            os.mkdir(gal_dir+band+'/Backsub_Temp')
            os.mkdir(gal_dir+band+'/SWarp_Temp')

            # Create Montage FITS header
            location_string = str(ra)+' '+str(dec)
            pix_size = bands_dict_req[band]['pix_size']
            montage_wrapper.commands.mHdr(location_string, width, gal_dir+band+'/'+str(name)+'_HDR', pix_size=pix_size)

            # Use Montage wrapper to reproject all fits files to common projection, skipping if none acually overlap
            print 'Performing reporjections for '+name+'_Spitzer_'+band+' maps'
            location_string = str(ra)+' '+str(dec)
            target_files = []
            proj_fail = 0
            [ target_files.append(target_file) for target_file in os.listdir(gal_dir+band) if '.fits' in target_file ]
            for target_file in target_files:
                try:
                    montage_wrapper.wrappers.reproject(os.path.join(gal_dir+band,target_file), os.path.join(gal_dir+band,target_file), header=gal_dir+band+'/'+str(name)+'_HDR', exact_size=True)
                except:
                    os.remove(os.path.join(gal_dir+band,target_file))
                    proj_fail += 1
            if proj_fail==len(target_files):
                print 'No Spitzer coverage for '+name+' at '+band
                continue

            # Loop over error maps and copy
            for listfile in os.listdir(gal_dir+band):
                if '_munc.fits' in listfile:
                    shutil.copy2(gal_dir+band+'/'+listfile, gal_dir+'Errors/'+band)

                    # Convert error maps to weight maps
                    unc_fitsdata = astropy.io.fits.open(gal_dir+band+'/'+listfile)
                    unc_image = unc_fitsdata[0].data
                    unc_header = unc_fitsdata[0].header
                    unc_fitsdata.close()
                    unc_image = unc_image**-1.0
                    unc_hdu = astropy.io.fits.PrimaryHDU(data=unc_image, header=unc_header)
                    unc_hdulist = astropy.io.fits.HDUList([unc_hdu])
                    unc_hdulist.writeto(gal_dir+band+'/SWarp_Temp/'+listfile.replace('_munc.fits','_maic.wgt.fits'), clobber=True)

                    # Delete old uncertainty map
                    os.remove(gal_dir+band+'/'+listfile)

            # If only one image file, proceed straight to co-adding; otherwise, commence background-matching
            mosaic_count = 0
            for listfile in os.listdir(gal_dir+band):
                if '.fits' in listfile:
                    mosaic_count += 1
            if mosaic_count==1:
                for listfile in os.listdir(gal_dir+band):
                    if '.fits' in listfile:
                        shutil.move(listfile, gal_dir+band+'/SWarp_Temp')
            if mosaic_count>1:

                # Use Montage wrapper to determine appropriate corrections for background matching
                print 'Determining background corrections for '+name+'_Spitzer_'+band+' maps'
                montage_wrapper.commands.mImgtbl(gal_dir+band,  gal_dir+band+'/'+band+'_Image_Metadata_Table.dat', corners=True)
                montage_wrapper.commands.mOverlaps(gal_dir+band+'/'+band+'_Image_Metadata_Table.dat', gal_dir+band+'/'+band+'_Image_Diffs_Table.dat')
                montage_wrapper.commands.mDiffExec(gal_dir+band+'/'+band+'_Image_Diffs_Table.dat', gal_dir+band+'/'+str(name)+'_HDR', gal_dir+band+'/Diffs_Temp', no_area=True)
                montage_wrapper.commands.mFitExec(gal_dir+band+'/'+band+'_Image_Diffs_Table.dat', gal_dir+band+'/'+band+'_Image_Fitting_Table.dat', gal_dir+band+'/Diffs_Temp')
                montage_wrapper.commands.mBgModel(gal_dir+band+'/'+band+'_Image_Metadata_Table.dat', gal_dir+band+'/'+band+'_Image_Fitting_Table.dat', gal_dir+band+'/'+band+'_Image_Corrections_Table.dat', level_only=True, n_iter=16384)

                # Apply background corrections using Montage subprocess, with timeout handling
                print 'Applying background corrections to '+name+'_Spitzer_'+band+' maps'
                mBgExec_fail_count = 0
                mBgExec_success = False
                mBgExec_uberfail = False
                while mBgExec_success==False:

                    # Attempt background-matching
                    mBgExec_sp = subprocess.Popen( ['/home/soft/montage/bin/mBgExec', '-n', gal_dir+band+'/'+band+'_Image_Metadata_Table.dat', gal_dir+band+'/'+band+'_Image_Corrections_Table.dat', gal_dir+band+'/SWarp_Temp' ], preexec_fn=os.setsid, stdout=subprocess.PIPE )
                    mBgExec_fail = False
                    seconds = 0
                    minutes_max = 45
                    while mBgExec_fail==False:
                        time.sleep(1)
                        mBgExec_stdout = mBgExec_sp.stdout.readline()
                        if mBgExec_sp.poll()==None:
                            seconds += 1
                        if 'Table has no data records' in mBgExec_stdout:
                            mBgExec_fail = True
                            mBgExec_fail_count += 1
                            break
                        if seconds>=(60*minutes_max):
                            mBgExec_fail = True
                            mBgExec_fail_count += 1
                            break
                        if mBgExec_sp.poll()!=None:
                            mBgExec_success = True
                            break

                    # Handle timeouts and other failures
                    if mBgExec_fail_count>0:
                        print 'Background matching with Montage has failed '+str(mBgExec_fail_count)+' time(s); reattempting'
                    if mBgExec_fail==True and mBgExec_success==False and mBgExec_fail_count>=3:
                        mBgExec_uberfail = True
                        print 'Background matching with Montage has failed 3 times; proceeding directly to co-additon'
                        try:
                            os.killpg( os.getpgid(mBgExec_sp.pid), 15 )
                        except:
                            'Background matching subprocess appears to have imploded; no task to kill'
                        for listfile in os.listdir(gal_dir+band):
                            if '_maic.fits' in listfile:
                                shutil.move(listfile, gal_dir+band+'/SWarp_Temp')
                        break



            # Sort out daft filename differences between image maps and error maps
            for gal_file in os.listdir(gal_dir+band+'/SWarp_Temp'):
                os.rename(gal_dir+band+'/SWarp_Temp/'+gal_file, gal_dir+band+'/SWarp_Temp/'+gal_file.replace('_'+gal_file.split('_')[-2:][0], '') )

            # Perform least-squares plane fitting to match MIPS image levels
            if instrument=='MIPS':
                ChrisFuncs.Coadd.LevelFITS(gal_dir+band+'/SWarp_Temp', 'maic.fits', convfile_dir=False)

            # Use SWarp to co-add images weighted by their error maps
            print 'Co-adding '+name+'_Spitzer_'+band+' maps'
            os.chdir(gal_dir+band+'/SWarp_Temp')
            os.system('swarp *_maic.fits -IMAGEOUT_NAME '+name+'_Spitzer_'+band+'_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -GAIN_KEYWORD DIESPIZERDIE -RESCALE_WEIGHTS N -SUBTRACT_BACK N -RESAMPLE N -VMEM_MAX 4095 -MEM_MAX 4096 -WEIGHT_TYPE MAP_WEIGHT -NTHREADS 4 -VERBOSE_TYPE QUIET')
            Spitzer_SWarp_NaN(name+'_Spitzer_'+band+'_SWarp.fits')

            # Re-project finalised image map using Montage
            montage_wrapper.wrappers.reproject(gal_dir+band+'/SWarp_Temp/'+name+'_Spitzer_'+band+'_SWarp.fits', out_dir+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'.fits', header=gal_dir+band+'/'+str(name)+'_HDR', exact_size=True)

            # Compress finalised image map
            os.chdir(out_dir)
            if os.path.exists(out_dir+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'.fits.gz'):
                os.remove(out_dir+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'.fits.gz')
            os.system('gzip '+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'.fits')
            print 'Completed Montaging and SWarping '+name+'_Spitzer_'+band+' image map'



            # Turn error maps into exposure time maps
            for listfile in os.listdir(gal_dir+'Errors/'+band):
                if '_munc.fits' in listfile:
                    unc_fitsdata = astropy.io.fits.open(gal_dir+'Errors/'+band+'/'+listfile)
                    unc_image = unc_fitsdata[0].data
                    unc_header = unc_fitsdata[0].header
                    unc_fitsdata.close()
                    unc_image = unc_image**-2.0
                    unc_hdu = astropy.io.fits.PrimaryHDU(data=unc_image, header=unc_header)
                    unc_hdulist = astropy.io.fits.HDUList([unc_hdu])
                    unc_hdulist.writeto(gal_dir+'Errors/'+band+'/'+listfile.replace('_munc.fits','_mexp.fits'), clobber=True)

            # Use Montage to add exposure time images
            print 'Co-adding '+name+'_Spitzer_'+band+' error maps'
            target_files = []
            [ target_files.append(dir_file) for dir_file in os.listdir(gal_dir+'Errors/'+band) if 'mexp.fits' in dir_file ]
            for i in range(0, len(target_files)):
                exp_fitsdata = astropy.io.fits.open(gal_dir+'Errors/'+band+'/'+target_files[i])
                exp_image = exp_fitsdata[0].data
                exp_header = exp_fitsdata[0].header
                exp_fitsdata.close(gal_dir+'Errors/'+band+'/'+target_files[i])
                if i==0:
                    add_image = np.zeros([ exp_image.shape[0], exp_image.shape[1] ])
                    add_header = exp_header.copy()
                exp_good = np.where( np.isnan(exp_image)==False )
                add_image[exp_good] += exp_image[exp_good]
            add_hdu = astropy.io.fits.PrimaryHDU(data=add_image, header=add_header)
            add_hdulist = astropy.io.fits.HDUList([add_hdu])
            add_hdulist.writeto(gal_dir+'Errors/'+band+'/'+name+'_Spitzer_'+band+'_Exp_Add.fits', clobber=True)

            # Re-project final exposure map using Montage
            montage_wrapper.wrappers.reproject(gal_dir+'Errors/'+band+'/'+name+'_Spitzer_'+band+'_Exp_Add.fits', gal_dir+'Errors/'+band+'/'+name+'_Spitzer_'+band+'_Exp.fits', header=gal_dir+band+'/'+str(name)+'_HDR', exact_size=True)

            # Convert final exposure time map into error map
            unc_fitsdata = astropy.io.fits.open(gal_dir+'Errors/'+band+'/'+name+'_Spitzer_'+band+'_Exp.fits')
            unc_image = unc_fitsdata[0].data
            unc_header = unc_fitsdata[0].header
            unc_fitsdata.close()
            unc_image[ np.where(unc_image<0) ] = np.NaN
            unc_image = unc_image**-0.5
            unc_image[ np.where(unc_image==np.inf) ] = np.NaN
            unc_hdu = astropy.io.fits.PrimaryHDU(data=unc_image, header=unc_header)
            unc_hdulist = astropy.io.fits.HDUList([unc_hdu])
            unc_hdulist.writeto(out_dir+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'_Error.fits', clobber=True)

            # Compress finalised exposure time map
            os.chdir(out_dir)
            if os.path.exists(out_dir+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'_Error.fits.gz'):
                os.remove(out_dir+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'_Error.fits.gz')
            os.system('gzip '+name+'_Spitzer_'+bands_dict_req[band]['band_long']+'_Error.fits')
            print 'Completed Montaging '+name+'_Spitzer_'+band+' error map'




        # Record that processing of souce has been compelted
        already_processed_file = open(already_processed_path, 'a')
        already_processed_file.write(name+'\n')
        already_processed_file.close()

        # Clean memory, and return timings
        shutil.rmtree(gal_dir)
        time_list.append( time.time() )
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        time_file = open( os.path.join('/'.join(in_dir.split('/')[:-2]),'Estimated_Completion_Time.txt'), 'w')
        time_file.write(time_est)
        time_file.close()
        print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'
