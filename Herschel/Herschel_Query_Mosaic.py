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



# Define function to wget and extact Herschel files
def Herschel_wget(data_url, data_filename):
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
            time.sleep(0.1)
            success = False
        os.system('gzip -d '+data_filename)



# Define function to replace null pixels in SWarp outputs with NaNs
def Herschel_SWarp_NaN(target):
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
if __name__  ==  "__main__":

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Herschel/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Herschel/Mosaics/'

    # Read in source catalogue
    ness_cat = np.genfromtxt(os.path.join(dropbox,'Work/Tables/NESS/NESS_Sample.csv'), delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']

    # State desired map width (in degrees)
    width = 0.5

    # Read in list of already-Montaged sources
    already_file = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/Herschel/Herschel_Already_Processed_List.dat'
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
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # State band information
    bands_dict = {'70':{'band':'70','instrument':'PACS','wavelength':'70um','filter':'PHOTBLUE','pix_size':2,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTBLUE','hdr_err_ext_name':'stDev'},
                  '100':{'band':'100','instrument':'PACS','wavelength':'100um','filter':'PHOTGREEN','pix_size':3,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTGREEN','hdr_err_ext_name':'stDev'},
                  '160':{'band':'160','instrument':'PACS','wavelength':'160um','filter':'PHOTRED','pix_size':4,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTRED','hdr_err_ext_name':'stDev'},
                  '250':{'band':'250','instrument':'SPIRE','wavelength':'250um','filter':'PSW','pix_size':6,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PSW','hdr_err_ext_name':'error'},
                  '350':{'band':'350','instrument':'SPIRE','wavelength':'350um','filter':'PMW','pix_size':8,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PMW','hdr_err_ext_name':'error'},
                  '500':{'band':'500','instrument':'SPIRE','wavelength':'500um','filter':'PLW','pix_size':12,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PLW','hdr_err_ext_name':'error'}}

    # State map mode prefixes we care about
    req_obs_modes = ['SpirePhotoLargeScan','SpirePhotoSmallScan','PacsPhoto','SpirePacsParallel']

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
        print 'Processing source '+name

        # Check if source is in list of already-montaged sources
        if name in already_processed:
            print name+' already processed'
            continue

        # Create data processing dirctory (deleting any prior), and set appropriate Python (ie, Montage) working directory
        gal_dir = in_dir+str(name)+'/'
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)
        os.makedirs(os.path.join(gal_dir,'Raw'))
        os.chdir(os.path.join(gal_dir,'Raw'))

        # Create band-specific directories
        for band in bands_dict.keys():
            os.makedirs(os.path.join(gal_dir,'Raw',band))



        # Perform query, with error handling
        print 'Querying HSA'
        query_success = False
        query_fail_count = 0
        while query_success == False:
            if query_fail_count>=10:
                break
            try:
                query_url = 'http://archives.esac.esa.int/hsa/aio/jsp/siap.jsp?POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS'
                query_filename = os.path.join(in_dir,name,str(name)+'.vot')
                if os.path.exists(query_filename):
                    os.remove(query_filename)
                wget.download(query_url, out=query_filename)
                query_success = True
            except:
                print 'HSA query failed; reattempting'
                query_fail_count += 1
        if not os.path.exists(query_filename):
            query_success=False

        # Read query result VOTable
        query_output = astropy.io.votable.parse_single_table(query_filename)
        query_table = query_output.array

        # If query finds no matches, continue to next target
        if len(query_table) == 0:
            print 'No Herschel data for '+name
            alrady_processed_file = open(already_file, 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()
            shutil.rmtree( os.path.join( in_dir, name ) )
            time_list.append(time.time())
            continue

        # Record which urls correspond to data in the desired modes (dealing with awkwardness for if there is only 1 entry, or silly massive files)
        hsa_urls = []
        if query_table.size == 1:
            if query_table['OBS_MODE'] in req_obs_modes:
                hsa_urls.append(query_table['DATA_ACCESS'])
        else:
            for j in range(0, query_table.size):
                if query_table['OBS_MODE'][j] in req_obs_modes:
                    hsa_urls.append(query_table['DATA_LINK'][j])

        # In parallel, download and extract files
        os.chdir(os.path.join(gal_dir,'Raw'))
        dl_pool = mp.Pool(processes=20)
        for j in range(0, len(hsa_urls)):
            data_url = hsa_urls[j]
            data_filename = os.path.join(gal_dir,'Raw',name+'_'+str(j)+'_HSA.fits')
            dl_pool.apply_async( Herschel_wget, args=(data_url, data_filename,) )
            #Herschel_wget(data_url, data_filename)
        dl_pool.close()
        dl_pool.join()

        # Loop over bands, and downloaded files (skipping folders), to move files to band-specific directories
        for band in bands_dict.keys():
            for listfile in os.listdir(os.path.join(gal_dir,'Raw')):
                if '.fits' not in listfile:
                    continue
                list_hdr = astropy.io.fits.getheader( os.path.join(gal_dir,'Raw',listfile), ext=0 )
                if list_hdr['INSTRUME'] == bands_dict[band]['instrument']:
                    if list_hdr[bands_dict[band]['hdr_inst_card_kwrd']] == bands_dict[band]['hdr_inst_card_entry']:
                        shutil.copy2(os.path.join(gal_dir,'Raw',listfile), os.path.join(gal_dir,'Raw',band))
                        os.remove(os.path.join(gal_dir,'Raw',listfile))

        # Loop over PACS bands and files to delete dud PACS calibration(?) maps
        for band in bands_dict.keys():
            if bands_dict[band]['instrument'] == 'PACS':
                for listfile in os.listdir(os.path.join(gal_dir,'Raw',band)):
                    if astropy.io.fits.getheader( os.path.join(gal_dir,'Raw',band,listfile), ext=0 )['OBSERVER'][-4:].lower() == 'pacs':
                        os.remove(os.path.join(gal_dir,'Raw',band,listfile))

        # Loop over each band's files, to save image map to separate FITS files
        for band in bands_dict.keys():
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band)):
                img_map, img_header = astropy.io.fits.getdata( os.path.join(gal_dir,'Raw',band,listfile), header=True, extname='image')

                # Record which image pixels are zeros, and convert to NaNs
                where_zero = np.where(img_map==0)
                img_map[where_zero] = np.NaN
                astropy.io.fits.writeto(os.path.join(gal_dir,'Raw',band,listfile.replace('.fits','_Img.fits')), img_map, header=img_header)

                # Now save coverage and error maps to separate files, with zeros similarly converted to NaNs
                cov_map, cov_header = astropy.io.fits.getdata( os.path.join(gal_dir,'Raw',band,listfile), header=True, extname='coverage')
                cov_map[where_zero] = np.NaN
                astropy.io.fits.writeto(os.path.join(gal_dir,'Raw',band,listfile.replace('.fits','_Cov.fits')), cov_map, header=cov_header)
                err_map, err_header = astropy.io.fits.getdata( os.path.join(gal_dir,'Raw',band,listfile), header=True, extname=bands_dict[band]['hdr_err_ext_name'])
                err_map[where_zero] = np.NaN
                astropy.io.fits.writeto(os.path.join(gal_dir,'Raw',band,listfile.replace('.fits','_Error.fits')), img_map, header=err_header)

        # Loop over each band for coaddition
        for band in bands_dict.keys():
            if not os.path.exists(os.path.join(gal_dir,'Raw',band)):
                continue
            if len(os.path.join(gal_dir,'Raw',band)) == 0:
                continue
            print 'Commencing Montaging and SWarping of '+name+'_Herschel_'+band

            # Create processing directories
            os.chdir(os.path.join(gal_dir,'Raw',band))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Img_Maps'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Cov_Maps'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Err_Maps'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Exp_Maps'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Wgt_Temp'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Pff_Temp'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'Backsub_Temp'))
            os.mkdir(os.path.join(gal_dir,'Raw',band,'SWarp_Temp'))

            # Create Montage FITS header
            location_string = str(ra)+' '+str(dec)
            pix_size = bands_dict[band]['pix_size']
            montage_wrapper.commands.mHdr(location_string, width, os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), pix_size=pix_size)

            # Use Montage wrapper to reproject all fits files to common projection, skipping if none acually overlap
            print 'Performing reporjections for '+name+'_Herschel_'+band+' maps'
            target_files = []
            proj_fail = 0
            [ target_files.append(target_file) for target_file in os.listdir(os.path.join(gal_dir,'Raw',band)) if '.fits' in target_file ]
            for target_file in target_files:
                try:
                    montage_wrapper.wrappers.reproject(os.path.join(os.path.join(gal_dir,'Raw',band,target_file)),
                                                       os.path.join(os.path.join(gal_dir,'Raw',band,target_file)),
                                                       header=os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'),
                                                       exact_size=True)
                except:
                    os.remove(os.path.join(os.path.join(gal_dir,'Raw',band,target_file)))
                    proj_fail += 1
            if proj_fail == len(target_files):
                print 'No Herschel coverage for '+name+' at '+band
                continue

            # Move reprojcted maps to relevant locations
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band)):
                if '_Img.fits' in os.path.join(gal_dir,'Raw',band,listfile):
                    shutil.move(os.path.join(gal_dir,'Raw',band,listfile), os.path.join(gal_dir,'Raw',band,'Img_Maps'))
                elif '_Cov.fits' in os.path.join(gal_dir,'Raw',band,listfile):
                    shutil.move(os.path.join(gal_dir,'Raw',band,listfile), os.path.join(gal_dir,'Raw',band,'Cov_Maps'))
                elif '_Error.fits' in os.path.join(gal_dir,'Raw',band,listfile):
                    shutil.move(os.path.join(gal_dir,'Raw',band,listfile), os.path.join(gal_dir,'Raw',band,'Err_Maps'))



            # If only one image file, proceed straight to co-adding; otherwise, commence background-matching
            mosaic_count = 0
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'Img_Maps')):
                if '_Img.fits' in listfile:
                    mosaic_count += 1
            if mosaic_count == 1:
                for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'Img_Maps')):
                    if '.fits' in listfile:
                        shutil.move(os.path.join(gal_dir,'Raw',band,'Img_Maps',listfile) , os.path.join(gal_dir,'Raw',band,'SWarp_Temp'))
            if mosaic_count>1:

                # Use Montage wrapper to determine appropriate corrections for background matching
                print 'Determining background corrections for '+name+'_Herschel_'+band+' maps'
                os.chdir(os.path.join(gal_dir,'Raw',band,'Img_Maps'))
                montage_wrapper.commands.mImgtbl( os.path.join(gal_dir,'Raw',band,'Img_Maps'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), corners=True )
                montage_wrapper.commands.mOverlaps( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Diffs_Table.dat') )
                montage_wrapper.commands.mDiffExec( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Diffs_Table.dat'), os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), os.path.join(gal_dir,'Raw',band,'Pff_Temp'), no_area=True, proj_dir=os.path.join(gal_dir,'Raw',band,'Img_Maps'))
                montage_wrapper.commands.mFitExec( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Diffs_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Fitting_Table.dat'), os.path.join(gal_dir,'Raw',band,'Pff_Temp') )
                montage_wrapper.commands.mBgModel( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Fitting_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Corrections_Table.dat'), level_only=True, n_iter=16384)

                # Apply background corrections using Montage subprocess, with timeout handling
                print 'Applying background corrections to '+name+'_Herschel_'+band+' maps'
                mBgExec_fail_count = 0
                mBgExec_success = False
                mBgExec_uberfail = False
                while mBgExec_success == False:

                    # Attempt background-matching
                    mBgExec_sp = subprocess.Popen( ['/home/soft/montage/bin/mBgExec', '-n', '-p', os.path.join(gal_dir,'Raw',band,'Img_Maps'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Corrections_Table.dat'), os.path.join(gal_dir,'Raw',band,'SWarp_Temp') ], preexec_fn=os.setsid, stdout=subprocess.PIPE )
                    mBgExec_fail = False
                    seconds = 0
                    minutes_max = 45
                    while mBgExec_fail == False:
                        time.sleep(1)
                        mBgExec_stdout = mBgExec_sp.stdout.readline()
                        if mBgExec_sp.poll() == None:
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
                    if mBgExec_fail_count>5:
                        print 'Background matching with Montage has failed '+str(mBgExec_fail_count)+' time(s); reattempting'
                    if mBgExec_fail == True and mBgExec_success == False and mBgExec_fail_count>=3:
                        mBgExec_uberfail = True
                        print 'Background matching with Montage has failed 5 times; proceeding directly to co-additon'
                        try:
                            os.killpg( os.getpgid(mBgExec_sp.pid), 15 )
                        except:
                            'Background matching subprocess appears to have imploded; no task to kill'
                        break
            if mBgExec_uberfail:
                for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'Img_Maps')):
                    if '_HSA_Img.fits' in listfile:
                        shutil.move(listfile, os.path.join(gal_dir,'Raw',band,'SWarp_Temp'))



            # Create weight maps, and copy to SWarp directory
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'Cov_Maps')):
                if '.fits' in listfile:
                    shutil.copy2( os.path.join(gal_dir,'Raw',band,'Cov_Maps',listfile), os.path.join(gal_dir,'Raw',band,'SWarp_Temp') )
                    wgt_image, wgt_header = astropy.io.fits.getdata( os.path.join(gal_dir,'Raw',band,'Cov_Maps',listfile), header=True )
                    wgt_image = wgt_image**0.5
                    astropy.io.fits.writeto( os.path.join(gal_dir,'Raw',band,'SWarp_Temp',listfile.replace('_Cov.fits','_Wgt.fits')), wgt_image, header=wgt_header )

            # Sort out daft filename differences between image maps and error maps
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'SWarp_Temp')):
                os.rename( os.path.join(gal_dir,'Raw',band,'SWarp_Temp',listfile), os.path.join(gal_dir,'Raw',band,'SWarp_Temp',listfile.replace('_Img.fits','.fits')) )
            """
            # Perform least-squares plane fitting to match image levels
            ChrisFuncs.Coadd.LevelFITS(os.path.join(gal_dir,'Raw',band,'SWarp_Temp'), 'Img.fits', convfile_dir=False)
            """
            # Use SWarp to co-add images weighted by their coverage maps
            print 'Co-adding '+name+'_Herschel_'+band+' maps'
            os.chdir(os.path.join(gal_dir,'Raw',band,'SWarp_Temp'))
            os.system('swarp *HSA.fits -IMAGEOUT_NAME '+name+'_Herschel_'+band+'_SWarp.fits -WEIGHT_SUFFIX _Wgt.fits -WEIGHT_TYPE MAP_RMS -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -GAIN_KEYWORD DIESPIZERDIE -RESCALE_WEIGHTS N -SUBTRACT_BACK N -RESAMPLE N -VMEM_MAX 4095 -MEM_MAX 4096 -WEIGHT_TYPE MAP_WEIGHT -NTHREADS 4 -VERBOSE_TYPE QUIET')
            Herschel_SWarp_NaN(name+'_Herschel_'+band+'_SWarp.fits')

            # Check that the final maps provides actual coverage of the point in question
            coadd_image, coadd_header = astropy.io.fits.getdata(os.path.join(gal_dir,'Raw',band,'SWarp_Temp',name+'_Herschel_'+band+'_SWarp.fits'), header=True)
            coadd_wcs = astropy.wcs.WCS(coadd_header)
            coords_xy = np.round(coadd_wcs.all_world2pix(np.array([[ra, dec]]), 0)).astype(int)
            coord_i, coord_j = coords_xy[0,1], coords_xy[0,0]
            if np.isnan(np.nanmax(coadd_image[coord_i-2:coord_i+2+1, coord_j-2:coord_j+2+2])):
                print('No Herschel coverage for '+name+' at '+band)
                alrady_processed_file = open(already_file, 'a')
                alrady_processed_file.write(name+'\n')
                alrady_processed_file.close()
                shutil.rmtree( os.path.join( in_dir, name ) )
                time_list.append(time.time())
                continue

            # Re-project finalised image map using Montage
            montage_wrapper.wrappers.reproject(os.path.join(gal_dir,'Raw',band,'SWarp_Temp',name+'_Herschel_'+band+'_SWarp.fits'), os.path.join(out_dir,name+'_Herschel_'+band+'.fits'), header=os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), exact_size=True)

            # Compress finalised FITs file
            os.chdir(out_dir)
            os.system('gzip '+os.path.join(gal_dir,'Raw',band,name+'_Herschel_'+band+'.fits'))
            print 'Completed Montaging and SWarping '+name+'_Herschel_'+band+' image map'



            # Turn error maps into exposure time maps
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'Err_Maps')):
                if '_Error.fits' in listfile:
                    err_image, err_header = astropy.io.fits.getdata(os.path.join(gal_dir,'Raw',band,'Err_Maps',listfile), header=True)
                    err_image = err_image**-2.0
                    astropy.io.fits.writeto(os.path.join(gal_dir,'Raw',band,'Exp_Maps',listfile.replace('_Error.fits','_Exp.fits')), err_image, header=err_header)

            # Use Montage to add exposure time images
            print 'Co-adding '+name+'_Herschel_'+band+' error maps'
            target_files = []
            [ target_files.append(dir_file) for dir_file in os.listdir(os.path.join(gal_dir,'Raw',band,'Exp_Maps')) if '_Exp.fits' in dir_file ]
            for i in range(0, len(target_files)):
                exp_image, exp_header = astropy.io.fits.getdata(os.path.join(gal_dir,'Raw',band,'Exp_Maps',target_files[i]), header=True)
                if i == 0:
                    add_image = np.zeros([ exp_image.shape[0], exp_image.shape[1] ])
                    add_header = exp_header.copy()
                exp_good = np.where( np.isnan(exp_image) == False )
                add_image[exp_good] += exp_image[exp_good]
            add_hdu = astropy.io.fits.PrimaryHDU(data=add_image, header=add_header)
            add_hdulist = astropy.io.fits.HDUList([add_hdu])
            astropy.io.fits.writeto(os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp_Add.fits'), add_image, header=add_header, clobber=True)

            # Re-project final exposure map using Montage
            montage_wrapper.wrappers.reproject(os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp_Add.fits'), os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp.fits'), header=os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), exact_size=True)

            # Convert final exposure time map into error map
            err_image, err_header = astropy.io.fits.getdata(os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp.fits'), header=True)
            err_image[ np.where(err_image<0) ] = np.NaN
            err_image = err_image**-0.5
            err_image[ np.where(err_image == np.inf) ] = np.NaN
            astropy.io.fits.writeto(os.path.join(out_dir,name+'_Herschel_'+band+'_Error.fits'), err_image, header=err_header, clobber=True)



            # Compress finalised exposure time map
            os.chdir(out_dir)
            os.system('gzip '+os.path.join(gal_dir,'Raw',band,name+'_Herschel_'+band+'_Error.fits'))
            print 'Completed Montaging and SWarping '+name+'_Herschel_'+band+' image map'
            print 'Completed Montaging '+name+'_Herschel_'+band+' error map'



        # Record that processing of souce has been compelted
        already_processed_file = open(already_file, 'a')
        already_processed_file.write(name+'\n')
        already_processed_file.close()

        # Clean memory, and return timings
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        gc.collect()
        time_list.append( time.time() )
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        time_file = open( os.path.join('/'.join(in_dir.split('/')[:-2]),'Estimated_Completion_Time.txt'), 'w')
        time_file.write(time_est)
        time_file.close()
        print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'
