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
import astropy.io.votable
import astropy.coordinates
from astroquery.skyview import SkyView
import gc
#warnings.filterwarnings("ignore")
import wget
import pdb
import time
import ChrisFuncs





# Define a timeout handler
def Handler(signum, frame):
    raise Exception("Timout!")



# Define function to wget DSS tiles
def Planck_wget(tile_url, tile_filename):
    if os.path.exists(tile_filename):
        os.remove(tile_filename)
    success = False
    while success==False:
        try:
            wget.download(tile_url, out=tile_filename)
            success = True
        except:
            time.sleep(0.1)
            success = False



# Define function to query for, and retrieve, DSS data from NASA SkyView
def Planck_SkyView(name, ra, dec, width, band, bands_dict, out_dir):
    print 'Retrieving '+band+' data for '+name+' from NASA SkyView'
    position_string = str(ra)+' '+str(dec)
    query_success = None
    query_reattempt = False
    pix_size = bands_dict[band]['pix_size']

    # Perform query
    while query_success==None:
        try:
            query_filename = out_dir+name+'_Planck_'+bands_dict[band]['wavelength']+'.fits'
            query_url = SkyView.get_image_list(position_string, bands_dict[band]['band_name'], deedger='_skip_', pixels=str(int((width*3600.0)/pix_size)), radius=astropy.units.Quantity(width, unit='deg'))
            if len(query_url)!=1:
                pdb.set_trace()
            query_success = True
        except:
            query_success = False

    # Retrieve and verify data
    if query_success:
        print 'Retrieving identified Planck '+band+' data for '+name
        Planck_wget(str(query_url[0]), query_filename)
        try:
            astropy.io.fits.info(query_filename)
        except Exception as exception:
            if exception.message=='Empty or corrupt FITS file':
                pdb.set_trace()
            else:
                pdb.set_trace()

    # Report failure
    if not query_success:
        print 'No Planck '+band+' data for '+name
        pdb.set_trace()





# Commence main task
if __name__ == "__main__":

    # Define paths
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Planck/Raw/'

    # Read in source catalogue
    ness_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Test_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = ness_cat['name']

    # Read in list of already-Montaged sources
    already_path = '/home/sarumandata2/spx7cjc/NESS/Test_Sample/Planck/Planck_Already_Processed_List.dat'
    if not os.path.exists(already_path):
        open(already_path,'a')
    already_processed = np.genfromtxt(already_path, dtype=None).tolist()

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

    # State band information
    bands_dict = {'Planck030':{'band_name':'Planck 030','wavelength':'10600','pix_size':204.0,'units':'K'},
                  'Planck044':{'band_name':'Planck 044','wavelength':'6810','pix_size':204.0,'units':'K'},
                  'Planck070':{'band_name':'Planck 070','wavelength':'4260','pix_size':204.0,'units':'K'},
                  'Planck100':{'band_name':'Planck 100','wavelength':'3000','pix_size':204.0,'units':'K'},
                  'Planck143':{'band_name':'Planck 143','wavelength':'2100','pix_size':102.0,'units':'K'},
                  'Planck217':{'band_name':'Planck 217','wavelength':'1380','pix_size':102.0,'units':'K'},
                  'Planck353':{'band_name':'Planck 353','wavelength':'850','pix_size':102.0,'units':'K'},
                  'Planck545':{'band_name':'Planck 545','wavelength':'550','pix_size':102.0,'units':'MJy/sr'},
                  'Planck857':{'band_name':'Planck 857','wavelength':'350','pix_size':102.0,'units':'MJy/sr'}}

    # Record time taken
    time_list = [time.time()]

    # Loop over each source
    for i in range(0, ness_cat.shape[0]):#np.where(name_list=='PGC042068')[0]:
        name = name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        width = 2.0
        print 'Processing target '+name

        # In parallel, retrieve DSS data in each band from NASA SkyView
        pool = mp.Pool(processes=5)
        for band in bands_dict.keys():
            pool.apply_async( Planck_SkyView, args=(name, ra, dec, width, band, bands_dict, out_dir,) )
            #Planck_SkyView(name, ra, dec, width, band, bands_dict, out_dir)
        pool.close()
        pool.join()

        # Record that processing of souce has been compelted
        alrady_processed_file = open(already_path, 'a')
        alrady_processed_file.write(name+'\n')
        alrady_processed_file.close()

        # Clean memory, and return timings
        gc.collect()
        time_list.append(time.time())
        time_est = ChrisFuncs.TimeEst(time_list, len(name_list))
        print 'Estimated completion time: '+time_est

# Jubilate
print 'All done!'



"""
query_fits = np.array(['.fits' in url for url in query_files])
query_tiles = list(query_files[np.where(query_fits==True)])
query_bands_dict = list(query_bands_dict[np.where(query_fits==True)])
"""
