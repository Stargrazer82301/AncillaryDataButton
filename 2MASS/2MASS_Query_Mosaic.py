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
import astropy
import astropy.io.votable
import montage_wrapper.commands
import cStringIO
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

# Define function to wget 2MASS tiles
def wget_2MASS(tile_url, tile_filename, verbose=False):
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

# Define function to Montage together contents of folder
def Montage_2MASS(name, ra, dec, pix_width, map_width, band, input_dir, out_dir):
    print 'Montaging '+name+'_2MASS_'+band
    if os.path.exists( os.path.join(input_dir,'Montage_Temp') ):
        shutil.rmtree( os.path.join(input_dir,'Montage_Temp') )
    location_string = str(ra)+' '+str(dec)
    os.chdir(input_dir)
    montage_wrapper.commands.mHdr(location_string, map_width, os.path.join(input_dir,str(name)+'_HDR'), pix_size=pix_width)
    montage_wrapper.wrappers.mosaic(input_dir, os.path.join(input_dir,'Montage_Temp'), header=None, n_proc=1, background_match=True, combine='mean', exact_size=False, level_only=True, cleanup=True, background_n_iter=32000 )
    montage_wrapper.wrappers.reproject(os.path.join(input_dir,'Montage_Temp','mosaic.fits'), os.path.join(out_dir,name+'_2MASS_'+band+'.fits'), header=os.path.join(input_dir,str(name)+'_HDR'), exact_size=True, cleanup=True, silent_cleanup=True)
    #os.system('mExec -r '+input_dir+' -f '+os.path.join(input_dir,str(name)+'_HDR')+' -o '+os.path.join(out_dir,name+'_2MASS_'+band+'.fits'))
    shutil.rmtree( os.path.join(input_dir,'Montage_Temp') )



# Commence main task
if __name__ == "__main__":

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/2MASS/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/2MASS/Mosaics/'

    # Read in source catalogue
    NESS_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = NESS_cat['name']

    # State bands of interest
    bands = ['J','H','K']
    vega_to_AB = [-0.91, -1.39, -1.85]

    # Delete any pre-existing temporary files
    if os.path.exists(in_dir):
        shutil.rmtree(in_dir)
    os.mkdir(in_dir)

    # Read in list of already-processed sources, and identify sources not yet processed
    already_file = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/2MASS/2MASS_Already_Processed_List.dat'
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
    for i in range(0, NESS_cat.shape[0]):#np.where(name_list=='NGC3031')[0]:
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        width = 0.25
        print 'Processing target '+name

        # Create tile processing dirctories
        os.mkdir( os.path.join( in_dir, name ) )
        os.mkdir( os.path.join( in_dir, name, 'Raw' ) )

        # Perform query (removing pre-existing query file, if present)
        print 'Querying 2MASS server using SIA protocol'
        query_success = False
        while not query_success:
            try:
                query_url = 'http://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS'
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
            print 'No 2MASS data for '+name
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
            if 'fits' not in query_table[j]['download']:
                continue
            tile_url = query_table[j]['download']
            tile_filename = os.path.join( in_dir, name, 'Raw', query_table[j]['download'].split('/')[-1:][0] )
            #wget_2MASS( tile_url, os.path.join(in_dir,name,'Raw',tile_filename) )
            dl_pool.apply_async( wget_2MASS, args=( tile_url, os.path.join(in_dir,name,'Raw',tile_filename), ) )
        dl_pool.close()
        dl_pool.join()

        # Create subfolders for each band
        for band in bands:
            os.mkdir( os.path.join( in_dir, name, band ) )

        # Copy files to relevant subfolders
        for j in range(0, len(query_table)):
            for raw_file in os.listdir( os.path.join( in_dir, name, 'Raw') ):
                if raw_file==query_table[j]['download'].split('/')[-1:][0]:
                    raw_band = query_table[j]['band']
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

            # Loop over maps, checking they are valid
            for raw_fits in os.listdir( os.path.join( in_dir, name, band) ):
                try:
                    raw_image, raw_header = astropy.io.fits.getdata( os.path.join( in_dir, name, band, raw_fits ), header=True )
                except Exception, exception_msg:
                    print exception_msg
                    print 'Invalid file: '+raw_fits
                    continue

                # Record pixel size
                raw_wcs = astropy.wcs.WCS(raw_header)
                raw_pix_width = np.mean( np.abs( raw_wcs.wcs.cdelt ) )
                if raw_pix_width<pix_min[b]:
                    pix_min[b] = raw_pix_width

                # Set all pixels to have units of Vega mag, using the header zero-point magnitude
                out_image = raw_header['MAGZP'] - ( 2.5 * np.log10(raw_image) )

                # Convert pixels value from Vega to AB magnitudes
                out_image = out_image - vega_to_AB[b]

                # Convert pixel units to Jy
                out_image = ChrisFuncs.ABMagsToJy(out_image)

                # Convert pixel units to Jy/sqdeg (switching to surface brightness here becuase different 2MASS surveys use different pixel sizes)
                pix_per_sqdeg = raw_pix_width**-2.0
                out_image *= pix_per_sqdeg

                # Make map zero mean
                out_image -= np.nanmean(out_image)

                # Save new image
                astropy.io.fits.writeto( os.path.join( in_dir, name, band, raw_fits ), out_image, header=raw_header, clobber=True )

        # In parallel, Montage together each band's input tiles, whilst dealing with timeouts
        out_pix_width = 1.0
        complete = False
        while not complete:
            try:
                signal.alarm(3600)
                pool = mp.Pool(processes=5)
                for b in range(0, len(bands)):
                    if bands_data[b]==False:
                        continue
                    band = bands[b]
                    pool.apply_async( Montage_2MASS, args=(name, ra, dec, out_pix_width, width, band, os.path.join( in_dir, name, band ), out_dir,) )
                    #Montage_2MASS(name, ra, dec, 0.4, width, band, os.path.join( in_dir, name, band ), out_dir)#3600.0*pix_min[b]
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
            if os.path.exists( os.path.join(out_dir,name+'_2MASS_'+band+'.fits') ):
                uncorr_image, uncorr_header = astropy.io.fits.getdata( os.path.join(out_dir,name+'_2MASS_'+band+'.fits'), header=True )
                uncorr_wcs = astropy.wcs.WCS(uncorr_header)
                uncorr_pix_width = np.mean( np.abs( uncorr_wcs.wcs.cdelt ) )
                sqdeg_per_pix = uncorr_pix_width**2.0
                corr_image = uncorr_image * sqdeg_per_pix
                astropy.io.fits.writeto( os.path.join(out_dir,name+'_2MASS_'+band+'.fits'), corr_image, header=uncorr_header, clobber=True )

        # Record that processing of souce has been compelted
        alrady_processed_file = open(already_file, 'a')
        alrady_processed_file.write(name+'\n')
        alrady_processed_file.close()

        # Clean memory, and return timings
        shutil.rmtree( os.path.join( in_dir, name ) )
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        time_file = open( os.path.join(in_dir,'Estimated_Completion_Time.txt'), 'w')
        time_file.write(time_est)
        time_file.close()
        print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'





"""
query_fits = np.array(['.fits' in url for url in query_files])
query_tiles = list(query_files[np.where(query_fits==True)])
query_bands = list(query_bands[np.where(query_fits==True)])
"""
