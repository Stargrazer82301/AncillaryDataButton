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
import subprocess
import astropy.io.fits
import astropy.wcs
import astropy.io.votable
import aplpy
import wget
import ChrisFuncs
plt.ioff()



# Define main function
def Run(ra, dec, width, name=None, out_dir=None, temp_dir=None, replace=False, flux=True, thumbnails=False,
        montage_path=None, swarp_path=None):
    """
    Function to generate standardised cutouts of Herschel observations.

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
                A string giving the path to be used as a temporary working directory by Herschel_Button. If not provided,
                a temporary directory will be created inside the output directory.
        replace: bool, optional
                If False, Herschel_Button will search the output directory for any pre-existing output FITS files from
                previous runs of the function, and will not bother repeat creating these maps (making it easy to resume
                processing a large number of targets from an interruption. If True, Herschel_Button will produce maps for
                all input targets, regardless of whether maps for these targets already exist in the output directory.
        flux: bool, optional
                If True, output maps will be in flux density units of Jy/pix. If false, output maps will be in surface
                brightness units of MJy/sr.
        thumbnails: bool, optional
                If True, JPG thumbnail images of the generated maps will also be proced and placed in out_dir.
        montage_path: str, optional
                Path to directory that contains the Montage commands (mProject, etc); useful if this directory is not in $PATH
        swarp_path: str: optional
                Path to directory that contains the SWarp command; useful if this directory is not in $PATH
    """


    # Handle Montage and SWarp paths, if kwargs provided
    if montage_path != None:
        os.environ['PATH'] += ':'+montage_path
    if swarp_path != None:
        os.environ['PATH'] += ':'+swarp_path
    import montage_wrapper

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

    # State band information
    bands_dict = {'70':{'band':'70','instrument':'PACS','wavelength':'70um','filter':'PHOTBLUE','pix_size':2,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTBLUE','hdr_err_ext_name':'stDev'},
                  '100':{'band':'100','instrument':'PACS','wavelength':'100um','filter':'PHOTGREEN','pix_size':3,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTGREEN','hdr_err_ext_name':'stDev'},
                  '160':{'band':'160','instrument':'PACS','wavelength':'160um','filter':'PHOTRED','pix_size':4,'hdr_inst_card_kwrd':'CAMERA','hdr_inst_card_entry':'PHOTRED','hdr_err_ext_name':'stDev'},
                  '250':{'band':'250','instrument':'SPIRE','wavelength':'250um','filter':'PSW','pix_size':6,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PSW','hdr_err_ext_name':'error'},
                  '350':{'band':'350','instrument':'SPIRE','wavelength':'350um','filter':'PMW','pix_size':8,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PMW','hdr_err_ext_name':'error'},
                  '500':{'band':'500','instrument':'SPIRE','wavelength':'500um','filter':'PLW','pix_size':12,'hdr_inst_card_kwrd':'DETECTOR','hdr_inst_card_entry':'PLW','hdr_err_ext_name':'error'}}

    # State map mode prefixes we care about
    req_obs_modes = ['SpirePhotoLargeScan','SpirePhotoSmallScan','PacsPhoto','SpirePacsParallel']



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
                if os.path.exists(os.path.join(out_dir,name+'_Herschel_'+bands_dict[band]['wavelength']+'.fits.gz')):
                    bands_done += 1

                # Also check for null files, indicated data not available for a givne band
                elif os.path.exists(os.path.join(out_dir,'.'+name+'_Herschel_'+bands_dict[band]['wavelength']+'.null')):
                    bands_done += 1

            # If this source has already been processed in all bands, skip it
            if bands_done == len(bands_dict.keys()):
                print('Herschel data for '+name+ ' already processed (if available); continuing to next target')
                time_list.append(time_list)
                continue
        print('Processing Herschel data for target '+name)

        # Create field processing dirctories (deleting any prior)
        gal_dir = os.path.join(temp_dir,str(name))+'/'
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)
        os.makedirs(os.path.join(gal_dir,'Raw'))
        os.chdir(os.path.join(gal_dir,'Raw'))

        # Create band-specific directories
        for band in bands_dict.keys():
            os.makedirs(os.path.join(gal_dir,'Raw',band))



        # Perform query, with error handling
        print('Querying HSA')
        query_success = False
        query_fail_count = 0
        while query_success == False:
            if query_fail_count>=10:
                break
            try:
                query_url = 'http://archives.esac.esa.int/hsa/aio/jsp/siap.jsp?POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS'
                query_filename = os.path.join(temp_dir,name,str(name)+'.vot')
                if os.path.exists(query_filename):
                    os.remove(query_filename)
                wget.download(query_url, out=query_filename)
                query_success = True
            except:
                print('HSA query failed; reattempting')
                query_fail_count += 1
        if not os.path.exists(query_filename):
            query_success=False

        # Read query result VOTable
        query_output = astropy.io.votable.parse_single_table(query_filename)
        query_table = query_output.array

        # Check if query returned any results; if not, create null file, and continue to next target
        if len(query_table)==0:
            print('No Herschel coverage for '+name+'; continuing to next target')
            os.system('touch '+os.path.join(temp_dir,'.'+name+'_Herschel_'+band+'.null'))
            continue

        # Record which urls correspond to data in the desired modes (dealing with awkwardness for if there is only 1 entry, or silly massive files)
        hsa_urls = []
        if query_table.size == 1:
            if query_table['OBS_MODE'] in req_obs_modes:
                hsa_urls.append(query_table['DATA_ACCESS'])
        else:
            for j in range(0, query_table.size):
                if query_table['OBS_MODE'][j].decode('utf-8') in req_obs_modes:
                    hsa_urls.append(query_table['DATA_LINK'][j].decode('utf-8'))

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
            print('Commencing Montaging and SWarping of '+name+'_Herschel_'+band)

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
            montage_wrapper.mHdr(location_string, width, os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), pix_size=pix_size)

            # Use Montage wrapper to reproject all fits files to common projection, skipping if none acually overlap
            print('Performing reprojections for '+name+'_Herschel_'+band+' maps')
            target_files = []
            proj_fail = 0
            [ target_files.append(target_file) for target_file in os.listdir(os.path.join(gal_dir,'Raw',band)) if '.fits' in target_file ]
            for target_file in target_files:
                try:
                    montage_wrapper.reproject(os.path.join(os.path.join(gal_dir,'Raw',band,target_file)),
                                              os.path.join(os.path.join(gal_dir,'Raw',band,target_file)),
                                              header=os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'),
                                              exact_size=True)
                except:
                    os.remove(os.path.join(os.path.join(gal_dir,'Raw',band,target_file)))
                    proj_fail += 1
            if proj_fail == len(target_files):
                print('No Herschel coverage for '+name+' at '+band)
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_Herschel_'+band+'.null'))
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
                print('Determining background corrections for '+name+'_Herschel_'+band+' maps')
                os.chdir(os.path.join(gal_dir,'Raw',band,'Img_Maps'))
                montage_wrapper.mImgtbl( os.path.join(gal_dir,'Raw',band,'Img_Maps'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), corners=True )
                montage_wrapper.mOverlaps( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Diffs_Table.dat') )
                montage_wrapper.mDiffExec( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Diffs_Table.dat'), os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), os.path.join(gal_dir,'Raw',band,'Pff_Temp'), no_area=True, proj_dir=os.path.join(gal_dir,'Raw',band,'Img_Maps'))
                montage_wrapper.mFitExec( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Diffs_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Fitting_Table.dat'), os.path.join(gal_dir,'Raw',band,'Pff_Temp') )
                montage_wrapper.mBgModel( os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Fitting_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Corrections_Table.dat'), level_only=True, n_iter=16384)

                # Apply background corrections using Montage subprocess, with timeout handling
                print('Applying background corrections to '+name+'_Herschel_'+band+' maps')
                mBgExec_fail_count = 0
                mBgExec_success = False
                mBgExec_uberfail = False
                while mBgExec_success == False:

                    # Attempt background-matching
                    mBgExec_sp = subprocess.Popen( ['mBgExec', '-n', '-p', os.path.join(gal_dir,'Raw',band,'Img_Maps'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Metadata_Table.dat'), os.path.join(gal_dir,'Raw',band,'Img_Maps',band+'_Image_Corrections_Table.dat'), os.path.join(gal_dir,'Raw',band,'SWarp_Temp') ], preexec_fn=os.setsid, stdout=subprocess.PIPE )
                    mBgExec_fail = False
                    seconds = 0
                    minutes_max = 45
                    while mBgExec_fail == False:
                        time.sleep(1)
                        mBgExec_stdout = mBgExec_sp.stdout.readline().decode()
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
                        print('Background matching with Montage has failed '+str(mBgExec_fail_count)+' time(s); reattempting')
                    if mBgExec_fail == True and mBgExec_success == False and mBgExec_fail_count>=3:
                        mBgExec_uberfail = True
                        print('Background matching with Montage has failed 5 times; proceeding directly to co-additon')
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
            print('Co-adding '+name+'_Herschel_'+band+' maps')
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
                os.system('touch '+os.path.join(temp_dir,'.'+name+'_Herschel_'+band+'.null'))
                continue

            # Re-project finalised image map using Montage
            montage_wrapper.reproject(os.path.join(gal_dir,'Raw',band,'SWarp_Temp',name+'_Herschel_'+band+'_SWarp.fits'), os.path.join(gal_dir,name+'_Herschel_'+band+'.fits'), header=os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), exact_size=True)

            # Compress finalised FITS file
            os.chdir(gal_dir)
            os.system('gzip '+os.path.join(gal_dir,name+'_Herschel_'+band+'.fits'))
            print('Completed processing '+name+'_Herschel_'+band+' image map')



            # Turn error maps into exposure time maps
            for listfile in os.listdir(os.path.join(gal_dir,'Raw',band,'Err_Maps')):
                if '_Error.fits' in listfile:
                    err_image, err_header = astropy.io.fits.getdata(os.path.join(gal_dir,'Raw',band,'Err_Maps',listfile), header=True)
                    err_image = err_image**-2.0
                    astropy.io.fits.writeto(os.path.join(gal_dir,'Raw',band,'Exp_Maps',listfile.replace('_Error.fits','_Exp.fits')), err_image, header=err_header)

            # Use Montage to add exposure time images
            print('Processing '+name+'_Herschel_'+band+' uncertainty map')
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
            montage_wrapper.reproject(os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp_Add.fits'), os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp.fits'), header=os.path.join(gal_dir,'Raw',band,str(name)+'.hdr'), exact_size=True)

            # Convert final exposure time map into error map
            err_image, err_header = astropy.io.fits.getdata(os.path.join(gal_dir,'Raw',band,'Exp_Maps',name+'_Herschel_'+band+'_Exp.fits'), header=True)
            err_image[ np.where(err_image<0) ] = np.NaN
            err_image = err_image**-0.5
            err_image[ np.where(err_image == np.inf) ] = np.NaN
            astropy.io.fits.writeto(os.path.join(gal_dir,name+'_Herschel_'+band+'_Error.fits'), err_image, header=err_header, clobber=True)

            # Compress finalised exposure time map
            os.chdir(out_dir)
            os.system('gzip '+os.path.join(gal_dir,name+'_Herschel_'+band+'_Error.fits'))
            print('Completed processing '+name+'_Herschel_'+band+' uncertainty map')



        # In parallel, generate final standardised maps for each band
        pool = mp.Pool(processes=9)
        for key in bands_dict.keys():
            band_dict = bands_dict[key]
            #pool.apply_async( Herschel_Generator, args=(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails,) )
            Herschel_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails)
        pool.close()
        pool.join()

        # Clean memory, and return timings (if more than one target being processed)
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        if len(name) > 1:
            print('Estimated time until Herschel data completed for all targets: '+time_est)

    # Tidy up (best as we can), and report completion
    gc.collect()
    try:
        shutil.rmtree(temp_dir)
    except:
        ChrisFuncs.RemoveCrawl(temp_dir)
        print('Unable to fully tidy up temporarny directory; probably due to NFS locks on network drive')
    print('All available Herschel imagery acquired for all targets')





# Define function to finalise Herschel image of a given source in a given band
def Herschel_Generator(name, ra, dec, temp_dir, out_dir, band_dict, flux, thumbnails):
    band = band_dict['band']
    wavelength = band_dict['wavelength']
    print('Generating final standardised map of Herschel '+band+' data for '+name)
    gal_dir = os.path.join(temp_dir,str(name))+'/'

    # If null file exists for this target in this band, copy it to final output directory
    if os.path.exists(os.path.join(temp_dir,'.'+name+'_Herschel_'+band+'.null')):
        shutil.copy(os.path.join(temp_dir,'.'+name+'_Herschel_'+band+'.null'),
                    os.path.join(out_dir,'.'+name+'_Herschel_'+band+'.null'))
    else:

        # Read in image map
        in_img, in_hdr = astropy.io.fits.getdata(os.path.join(gal_dir,name+'_Herschel_'+band+'.fits.gz'), header=True)
        in_wcs = astropy.wcs.WCS(in_hdr)
        in_pix_width_arcsec = 3600.0 * np.abs(np.max(in_wcs.pixel_scale_matrix))
        out_img = in_img.copy()

        # Read in error map
        in_unc = astropy.io.fits.getdata(gal_dir+name+'_Herschel_'+band+'_Error.fits.gz')
        out_unc = in_unc.copy()

        # Give default pixel units
        if band_dict['instrument'] == 'SPIRE':
            pix_unit = 'MJy/sr'
        if band_dict['instrument'] == 'PACS':
            pix_unit = 'Jy/pix'

        # If required, convert pixel units from MJy/sr to Jy/pix
        if (band_dict['instrument'] == 'SPIRE') and (flux == True):
            sr_per_sqarcsec = 2.3504E-11
            sr_per_pixels = sr_per_sqarcsec * in_pix_width_arcsec**2
            out_img *= 1E6 * sr_per_pixels
            out_unc *= 1E6 * sr_per_pixels
            pix_unit = 'Jy/pix'

        # If required, convert pixel units from MJy/sr to Jy/pix
        if (band_dict['instrument'] == 'PACS') and (flux == False):
            sr_per_sqarcsec = 2.3504E-11
            sr_per_pixels = sr_per_sqarcsec * in_pix_width_arcsec**2
            out_img *= 1E-6 * sr_per_pixels**-1
            out_unc *= 1E-6 * sr_per_pixels**-1
            pix_unit = 'MJy/sr'

        # Create standard header
        out_hdr = astropy.io.fits.Header()
        date = datetime.datetime.now().isoformat()

        # Populate standard header entries
        out_hdr.set('TARGET', name, 'Target source of this map')
        out_hdr.set('COORDSYS', 'IRCS', 'Coordinate reference frame for the RA and DEC')
        out_hdr.set('SIGUNIT', pix_unit, 'Unit of the map')
        out_hdr.set('TELESCOP', 'Herschel', 'Telescope that made this observation')
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

        # Create image-specific header and HDU
        out_hdr_img = out_hdr.copy()
        out_hdr_img.set('EXTNAME', 'image')
        image_out_hdu = astropy.io.fits.PrimaryHDU(data=out_img, header=out_hdr_img)

        # Create error-specific header and HDU
        out_hdr_unc = out_hdr.copy()
        out_hdr_unc.set('EXTNAME', 'error')
        error_out_hdu = astropy.io.fits.ImageHDU(data=out_unc, header=out_hdr_unc)

        # Create hdulist and save to file
        out_hdulist = astropy.io.fits.HDUList([image_out_hdu, error_out_hdu])
        out_hdulist.writeto(os.path.join(out_dir,name+'_Herschel_'+band+'.fits.gz'), overwrite=True)

        # Write output FITS file
        astropy.io.fits.writeto(os.path.join(out_dir,name+'_Herschel_'+band+'.fits.gz'), data=out_img, header=out_hdr, overwrite=True)

        # Make thumbnail image of cutout
        if thumbnails:
            try:
                thumb_out = aplpy.FITSFigure(os.path.join(out_dir,name+'_Herschel_'+band+'.fits'))
                thumb_out.show_colorscale(cmap='gist_heat', stretch='arcsinh')
                thumb_out.axis_labels.hide()
                thumb_out.tick_labels.hide()
                thumb_out.ticks.hide()
                thumb_out.show_markers(np.array([float(ra)]), np.array([float(dec)]), marker='+', s=500, lw=2.5, edgecolor='#01DF3A')
                thumb_out.save(os.path.join(out_dir,name+'_Herschel_'+band+'.jpg'), dpi=125)
                thumb_out.close()
            except:
                print('Failed making thumbnail for '+name)
                pdb.set_trace()

        # Clean memory before finishing
        gc.collect()





# Define function to wget and extact Herschel files
def Herschel_wget(data_url, data_filename):
    print('Acquiring '+data_url)
    if os.path.exists(data_filename):
        os.remove(data_filename)
    success = False
    while success == False:
        try:
            wget.download(data_url, out=data_filename)
            os.system('gzip -d '+data_filename)
            print('Successful acquisition of '+data_url)
            success = True
        except:
            print('Failure! Retrying acquistion of '+data_url)
            time.sleep(0.1)
            success = False





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





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")





#Run(23.4621, +30.6599, 1.0, name='M33', out_dir='/astro/dust_kg/cclark/Local_Dust/Raw_Obs/M33/Herschel/', montage_path='/Users/cclark/Soft/Montage/bin/', swarp_path='/usr/local/bin/')