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
import astropy.io.votable
import astropy.nddata.utils
import astropy.units
import astropy.coordinates
import montage_wrapper.commands
import shutil
import signal
import gc
import wget
import pdb
import time
import ChrisFuncs





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")

# Define function to wget VISTA tiles
def VISTA_wget(tile_url, tile_filename, verbose=True):
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

# Define function to Montage together contents of folder
def VISTA_Montage(name, ra, dec, pix_width, map_width, band, input_dir, out_dir):
    print 'Montaging '+name+'_VISTA_'+band
    if os.path.exists( os.path.join(input_dir,'Montage_Temp') ): shutil.rmtree( os.path.join(input_dir,'Montage_Temp') )
    if os.path.exists( os.path.join(input_dir,'Montage_Work') ): shutil.rmtree( os.path.join(input_dir,'Montage_Work') )
    location_string = str(ra)+' '+str(dec)
    os.chdir(input_dir)
    montage_wrapper.commands.mHdr(location_string, map_width, os.path.join(input_dir,str(name)+'_HDR'), pix_size=pix_width)
    montage_wrapper.wrappers.mosaic(input_dir, os.path.join(input_dir,'Montage_Temp'), header=None, background_match=True, level_only=False, combine='mean', exact_size=False, cleanup=True, background_n_iter=10000, work_dir=os.path.join(input_dir,'Montage_Work') )
    montage_wrapper.wrappers.reproject(os.path.join(input_dir,'Montage_Temp','mosaic.fits'), os.path.join(out_dir,name+'_VISTA_'+band+'.fits'), header=os.path.join(input_dir,str(name)+'_HDR'), exact_size=True, cleanup=True, silent_cleanup=True)
    #ChrisFuncs.FitsCutout( os.path.join(input_dir,'Montage_Temp','mosaic.fits'), ra, dec, 0.5*map_width*3600.0, reproj=False, variable=False, outfile=os.path.join(out_dir,name+'_VISTA_'+band+'.fits'), parallel=False )
    #os.system('mExec -r '+input_dir+' -f '+os.path.join(input_dir,str(name)+'_HDR')+' -o '+os.path.join(out_dir,name+'_VISTA_'+band+'.fits'))
    shutil.rmtree( os.path.join(input_dir,'Montage_Temp') )



# Commence main task
if __name__ == "__main__":

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/VISTA/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/VISTA/Mosaics/'

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']

    # State bands of interest
    bands = ['Z','Y','J','H','Ks']
    vega_to_AB = [-0.521, -0.618, -0.937, -1.384, -1.839]

    # Delete any pre-existing temporary files
    if os.path.exists(in_dir):
        shutil.rmtree(in_dir)
    os.mkdir(in_dir)

    # Read in list of already-processed sources, and identify sources not yet processed
    already_file = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/VISTA/VISTA_Already_Processed_List.dat'
    if not os.path.exists(already_file):
        open(already_file,'a')
    alrady_processed = np.genfromtxt(already_file, dtype=('S50')).tolist()
    remaining_list = []
    for i in range(0, name_list.shape[0]):
        already_done = 0
        name = name_list[i]
        if name not in alrady_processed:
            remaining_list.append(i)
    ness_cat = ness_cat[remaining_list]
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # Register signal function handler, for dealing with timeouts
    signal.signal(signal.SIGALRM, Handler)

    # Record time taken
    time_list = [time.time()]

    # Loop over each source
    for i in np.random.permutation(np.array(range(0, ness_cat.shape[0]))):
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        width = 0.25
        print 'Processing target '+name

        # Create tile processing dirctories
        os.mkdir( os.path.join( in_dir, name ) )
        os.mkdir( os.path.join( in_dir, name, 'Raw' ) )

        # Perform query (removing pre-existing query file, if present)
        print 'Querying VISTA server using SIA protocol'
        query_success = False
        while not query_success:
            try:
                query_url = 'http://wfaudata.roe.ac.uk/vista-siap/?POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS'
                query_filename = os.path.join(in_dir,name,str(name)+'.vot')
                if os.path.exists(query_filename):
                    os.remove(query_filename)
                wget.download(query_url, out=query_filename)
                query_success = True
            except:
                pdb.set_trace()
                query_success = False

        # Read query result VOTable
        query_table = astropy.table.Table.read(query_filename)

        # If query finds no matches, continue to next target
        if len(query_table)==0:
            print 'No VISTA data for '+name
            alrady_processed_file = open(already_file, 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()
            shutil.rmtree( os.path.join( in_dir, name ) )
            time_list.append(time.time())
            continue

        # Loop over files that are to be downloaded in parallel
        os.chdir(in_dir)
        dl_pool = mp.Pool(processes=20)
        for j in range(0, len(query_table)):

            # Skip those inpertinent tile products that start with 'e'
            if query_table[j]['Reference'].split('/')[-1:][0][0]=='e':
                continue

            # Alter URL so that full VISTA MEF file is retrieved, but still keep track of the extension of interest
            tile_extnum = query_table[j]['Reference'].split('http://surveys.roe.ac.uk/wsa/cgi-bin/siap.cgi?extNum=')[1].split('&file=')[0]
            tile_url = 'http://horus.roe.ac.uk/wsa/cgi-bin/fits_download.cgi?file='+query_table[j]['Reference'].split('file=')[1]+'&MFID='+query_table[j]['ObsId']
            tile_filename = os.path.join( in_dir, name, 'Raw', tile_extnum+'_'+tile_url.split('/')[-1:][0]+'.fits' )

            # Commence parallel download
            #VISTA_wget( tile_url, os.path.join(in_dir,name,'Raw',tile_filename), verbose=True )
            dl_pool.apply_async( VISTA_wget, args=( tile_url, os.path.join(in_dir,name,'Raw',tile_filename), ) )
        dl_pool.close()
        dl_pool.join()

        # Create subfolders for each band
        for band in bands:
            os.mkdir( os.path.join( in_dir, name, band ) )

        # Copy files to relevant subfolders
        for j in range(0, len(query_table)):
            for raw_file in os.listdir( os.path.join( in_dir, name, 'Raw' ) ):
                if '_'.join(raw_file.split('_')[1:]).split('&MFID')[0]==query_table[j]['Reference'].split('/')[-1:][0]:
                    raw_band = query_table[j]['Bandpass']
                    shutil.move( os.path.join( in_dir, name, 'Raw', raw_file), os.path.join( in_dir, name, raw_band) )
        shutil.rmtree( os.path.join( in_dir, name, 'Raw' ) )

        # Loop over bands
        pix_min = [np.inf]*len(bands)
        bands_data = [False]*len(bands)
        for b in range(0, len(bands)):
            band = bands[b]
            if len( os.listdir( os.path.join( in_dir, name, band) ) )>0:
                bands_data[b] = True
            else:
                continue

            # Loop over maps, first checking they are valid
            for raw_fits in os.listdir( os.path.join( in_dir, name, band) ):
                raw_exten = int(raw_fits.split('_')[0]) - 1
                try:
                    test_image, test_header = astropy.io.fits.getdata( os.path.join( in_dir, name, band, raw_fits ), header=True, ext=raw_exten )
                    print 'Preprocessing file: '+raw_fits
                except Exception, exception_msg:
                    print exception_msg
                    print 'Invalid file: '+raw_fits
                    continue

                # Load in appropriate file extensions, record pixel size, and retain original header for later use
                raw_fitsdata = astropy.io.fits.open( os.path.join( in_dir, name, band, raw_fits ) )
                raw_header = raw_fitsdata[raw_exten].header
                raw_image = raw_fitsdata[raw_exten].data
                raw_fitsdata.close()
                raw_wcs = astropy.wcs.WCS(raw_header)
                raw_pix_width = np.sqrt( np.max(np.abs(raw_wcs.pixel_scale_matrix))**2.0 + np.min(np.abs(raw_wcs.pixel_scale_matrix))**2.0 )
                if raw_pix_width<pix_min[b]:
                    pix_min[b] = raw_pix_width

                # Save data in more mundane FITS arrangement, and use this to produce cutout of relevant portion of map only, to speed processing
                try:
                    cutout_obj = astropy.nddata.utils.Cutout2D( raw_image, astropy.coordinates.SkyCoord(ra, dec, unit='deg'), 1.414*width*astropy.units.degree, wcs=raw_wcs, mode='trim', copy=True )
                except astropy.nddata.utils.NoOverlapError:
                    astropy.io.fits.writeto( os.path.join( in_dir, name, band, raw_fits ), np.zeros(raw_image.shape), header=raw_header, clobber=True )
                    continue
                in_image = cutout_obj.data
                in_wcs = cutout_obj.wcs
                in_header = in_wcs.to_header()

                # Use header information and Vega-to-AB conversions to work out true zero-point magnitude, as per prescription here: https://www.eso.org/sci/observing/phase3/data_releases/vvv_dr4.pdf
                header_exptime = raw_header['EXPTIME']
                header_airmass = 0.5 * ( raw_header['HIERARCH ESO TEL AIRM START'] + raw_header['HIERARCH ESO TEL AIRM END'] )
                header_extinct = raw_header['EXTINCT']
                header_mag_zpt = raw_header['MAGZPT']
                mag_zpt = header_mag_zpt - (2.5*np.log10(1/header_exptime)) - (header_extinct*(header_airmass-1.0)) - vega_to_AB[b]

                # Set all pixels to have units of AB mag, using the calculated zero-point magnitude
                out_image = mag_zpt - ( 2.5 * np.log10(in_image) )

                # Convert pixel units to Jy
                out_image = ChrisFuncs.ABMagsToJy(out_image)

                # Convert pixel units to Jy/sqdeg (switching to surface brightness here becuase different VISTA surveys use different pixel sizes)
                pix_per_sqdeg = raw_pix_width**-2.0
                out_image *= pix_per_sqdeg

                # Fix header WCS issue, and save new image
                in_header.set('CD1_1', in_header['PC1_1'] , 'Coordinate transformation matrix element')
                in_header.set('CD1_2', in_header['PC1_2'] , 'Coordinate transformation matrix element')
                in_header.set('CD2_1', in_header['PC2_1'] , 'Coordinate transformation matrix element')
                in_header.set('CD2_2', in_header['PC2_2'] , 'Coordinate transformation matrix element')
                in_header.remove('PC1_1')
                in_header.remove('PC1_2')
                in_header.remove('PC2_1')
                in_header.remove('PC2_2')
                astropy.io.fits.writeto( os.path.join( in_dir, name, band, raw_fits ), out_image, header=in_header, clobber=True )

        # In parallel, Montage together each band's input tiles, whilst dealing with timeouts
        complete = False
        fail_counter = 0
        while not complete:
            if fail_counter >= 5:
                continue
            try:
                signal.alarm(3600)
                pool = mp.Pool(processes=5)
                for b in range(0, len(bands)):
                    if bands_data[b]==False:
                        continue
                    band = bands[b]
                    pool.apply_async( VISTA_Montage, args=(name, ra, dec, 0.4, width, band, os.path.join( in_dir, name, band ), out_dir,) )
                    #VISTA_Montage(name, ra, dec, 0.4, width, band, os.path.join( in_dir, name, band ), out_dir)#3600.0*pix_min[b]
                pool.close()
                pool.join()
                complete = True
            except Exception, exception_msg:
                gc.collect()
                for band in bands:
                    input_dir = os.path.join(in_dir,name)+band+'/'
                    if os.path.exists(input_dir+'/Montage_Temp'):
                        shutil.rmtree(input_dir+'/Montage_Temp')
                print 'Mosaicing failure in '+band+'-band!'
                print exception_msg

        # Convert mosaiced map from units of Jy/sqdeg to Jy/pix
        for band in bands:
            if os.path.exists( os.path.join(out_dir,name+'_VISTA_'+band+'.fits') ):
                uncorr_image, uncorr_header = astropy.io.fits.getdata( os.path.join(out_dir,name+'_VISTA_'+band+'.fits'), header=True )
                uncorr_wcs = astropy.wcs.WCS(uncorr_header)
                uncorr_pix_width = np.mean( np.abs( uncorr_wcs.wcs.cdelt ) )
                sqdeg_per_pix = uncorr_pix_width**2.0
                corr_image = uncorr_image * sqdeg_per_pix
                astropy.io.fits.writeto( os.path.join(out_dir,name+'_VISTA_'+band+'.fits'), corr_image, header=uncorr_header, clobber=True )

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
        print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'
