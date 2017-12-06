# Identify location
import socket
location = socket.gethostname()
if location == 'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'Hobbitslayer':
    dropbox = 'C:\\Users\\spx7cjc\\Dropbox\\'
if location == 'saruman':
    dropbox = '/home/user/spx7cjc/Desktop/Herdata/Dropbox/'

# Import smorgasbord
import os
import sys
import multiprocessing as mp
import numpy as np
import astropy
import astropy.io.votable
import montage_wrapper.commands
import shutil
import signal
import gc
#warnings.filterwarnings("ignore")
import wget
import pdb
import time
import ChrisFuncs





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget WISE tiles
def WISE_wget(tile_url, tile_filename):
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



# Define function to Montage together error maps
def WISE_Montage(name, ra, dec, width, band, input_dir, out_dir):
    print 'Montaging '+name+'_WISE_'+band

    # Create header for reprojection
    location_string = str(ra)+' '+str(dec)
    montage_wrapper.commands.mHdr(location_string, width, os.path.join( input_dir, str(name)+'_HDR' ), pix_size=1.375)

    # Loop over error maps
    unc_files = os.listdir(input_dir)
    unc_first = True
    for i in range(0, len(unc_files)):
        unc_file = unc_files[i]
        if '.fits' not in unc_file:
            continue

        # Reproject error map, and read in
        montage_wrapper.wrappers.reproject( os.path.join( input_dir, unc_file ), os.path.join( input_dir, unc_file.replace('.fits.gz','_reproj.fits') ), header=os.path.join( input_dir, str(name)+'_HDR' ), exact_size=True, cleanup=True, silent_cleanup=True)
        unc_image, unc_header = astropy.io.fits.getdata( os.path.join( input_dir, unc_file.replace('.fits.gz','_reproj.fits') ), header=True )

        # Create array to co-add into
        if unc_first==True:
            wgt_coadd = np.zeros(unc_image.shape)
            unc_first = False

        # Convert error maps to exposure time maps, and co-add
        wgt_image = unc_image**-2.0
        wgt_bad = np.where( np.isnan(wgt_image)==True )
        wgt_image[ wgt_bad ] = 0.0
        wgt_coadd += wgt_image

    # Convert exposure time map back into error map, and save to file
    if unc_first==False:
        unc_coadd = wgt_coadd**-0.5
        astropy.io.fits.writeto( os.path.join( out_dir, name+'_WISE_'+band+'_Error.fits' ), unc_coadd, header=unc_header, clobber=True )





# Commence main task
if __name__ == "__main__":

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/WISE/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/WISE/Mosaics/'

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']

    # Read in list of already-Montaged sources
    already_processed_path = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/WISE/WISE_Already_Processed_List_Errors.dat'
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
    ness_cat = ness_cat[remaining_list]
    name_list = ness_cat['name']
    ra_list = ness_cat['ra']
    dec_list = ness_cat['dec']

    # State band infomration
    bands = ['W1','W2','W3','W4']

    # Register signal function handler, for dealing with timeouts
    signal.signal(signal.SIGALRM, Handler)

    # Record time taken
    time_list = [time.time()]

    # Loop over each source
    for i in range(0, ness_cat.shape[0]):#np.where(name_list=='NGC3031')[0]:
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        time_start = time.time()
        width = 0.25
        print 'Processing source '+name

        # Create tile processing dirctories (deleting any prior), and set appropriate Python (ie, Montage) working directory
        gal_dir = in_dir+str(name)+'/'
        if os.path.exists(gal_dir):
            shutil.rmtree(gal_dir)
        os.makedirs(gal_dir)
        os.makedirs(gal_dir+'/Raw')
        os.chdir(gal_dir+'/Raw')

        # Perform query (removing pre-existing query file, if present)
        print 'Querying IRSA server'
        query_success = False
        while not query_success:
            try:
                #query_url = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allsky/4band_p1bm_frm?POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS'
                query_url = 'http://irsa.ipac.caltech.edu/ibe/sia/wise/allwise/p3am_cdd?POS='+str(ra)+','+str(dec)+'&SIZE='+str(width)+'&INTERSECT=OVERLAPS'
                query_filename = gal_dir+str(name)+'.vot'
                if os.path.exists(query_filename):
                    os.remove(query_filename)
                wget.download(query_url, out=query_filename)
                query_success = True
            except:
                query_success = False

        # Read query result VOTable
        query_output = astropy.io.votable.parse_single_table(query_filename)
        query_table = query_output.array

        # Extract FITS file URLs, and their respective bands, from table
        query_urls = list( query_table['sia_url'] )
        query_fits_boolean = np.array([ '.fits' in url for url in query_urls ])
        query_urls = list( np.array(query_urls)[ np.where( query_fits_boolean==True ) ] )
        query_bands = list(query_table['sia_bp_id'])

        # Alter URLs so that they retrieve uncertainty maps instead
        query_urls = [ query_url.replace('-int-','-unc-').replace('.fits','.fits.gz') for query_url in query_urls ]

        # Delete VOTable
        os.remove(gal_dir+str(name)+'.vot')

        # In parallel, download files to holding directory
        dl_pool = mp.Pool(processes=20)
        for j in range(0, len(query_urls)):
            tile_url = query_urls[j]
            tile_filename = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/WISE/Temporary_Files/'+str(name)+'/Raw/'+str(name)+'_'+str(j)+'_'+str(query_bands[j])+'.fits.gz'
            #WISE_wget(tile_url, tile_filename)
            dl_pool.apply_async( WISE_wget, args=(tile_url, tile_filename,) )
        dl_pool.close()
        dl_pool.join()

        # Copy files to relevant folders
        for band in bands:
            if os.path.exists(gal_dir+band+'/'):
                shutil.rmtree(gal_dir+band+'/')
            os.makedirs(gal_dir+band+'/')
            raw_files = np.array(os.listdir(gal_dir+'Raw/'))
            for raw_file in raw_files:
                if band+'.fits' in raw_file:
                    shutil.move(gal_dir+'Raw/'+raw_file, gal_dir+band)
        shutil.rmtree(gal_dir+'Raw')

        # In parallel, Montage together each band's input tiles, whilst dealing with timeouts
        complete = False
        while not complete:
            try:
                signal.alarm(1800)
                pool = mp.Pool(processes=4)
                for band in bands:
                    WISE_Montage(name, ra, dec, width, band, os.path.join(gal_dir,band), out_dir)
                    #pool.apply_async( WISE_Montage, args=(name, ra, dec, width, band, os.path.join(gal_dir,band), out_dir,) )
                pool.close()
                pool.join()
                complete = True
            except Exception, exception_msg:
                gc.collect()
                for band in bands:
                    input_dir = gal_dir+band+'/'
                    if os.path.exists(input_dir+'/Montage_Temp'):
                        shutil.rmtree(input_dir+'/Montage_Temp')
                print exception_msg

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
