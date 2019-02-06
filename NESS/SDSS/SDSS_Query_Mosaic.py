# Identify location
import socket
location = socket.gethostname()
if location == 'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'Hobbitslayer':
    dropbox = 'C:\\Users\\spx7cjc\\Dropbox\\'
if location == 'saruman':
    dropbox = '/home/user/spx7cjc/Desktop/Herdata/Dropbox/'
if location == 'replicators.stsci.edu':
    dropbox = '/Users/cclark/Dropbox/'

# Import smorgasbord
import os
import sys
import multiprocessing as mp
import numpy as np
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
    raise Exception('Timout!')

# Define function to check if a list of SDSS DR12 URL correspond to a primary fields, and return only those that are
def SDSS_Primary_Check(urls, index):
    urls_pri = []
    for url in urls:
        run = url.split('/')[9].lstrip('0')
        camcol = url.split('/')[10].lstrip('0')
        field = url.split('/')[11].split('.')[0].split('-')[4].lstrip('0')
        check_string = run+' '+camcol+' '+field
        if check_string in index:
            urls_pri.append(url)
    return urls_pri

# Define function to consturct URL for a given SDSS field
def SDSS_URL(run, camcol, field, band):
    sdss_url = 'http://data.sdss3.org/sas/dr12/boss/photoObj/frames/301/'+run+'/'+camcol+'/frame-'+band+'-'+run.zfill(6)+'-'+camcol+'-'+field.zfill(4)+'.fits.bz2'
    return sdss_url
#
## Define function to wget SDSS fields; reject and re-acquire fields less than 1MB in size
#def SDSS_wget(tile_url, path):
#    fitsname = tile_url.split('/')[-1:][0]
#    print('Acquiring '+fitsname)
#    if os.path.exists(path+fitsname):
#        os.remove(path+fitsname)
#    else:
#        success = False
#        while success==False:
#            try:
#                wget.download(tile_url, out=path+fitsname)
#                filesize = os.stat(path+fitsname).st_size
#                if float(filesize)<1048576.0:
#                    raise NameError('File not large enough')
#                print('Successful acquisition of '+fitsname)
#                success = True
#            except:
#                print('Failure! Retrying acquistion of '+fitsname)
#                time.sleep(0.1)
#                success = False
#
## Define function to extract SDSS bz2 archives
#def SDSS_Extract(filepath):
#    print('Decompressing file '+str( filepath.split('/')[-1:][0] ))
#    os.system('bzip2 -d '+filepath)
#
#
#
## Define function to Montage together contents of folder
#def SDSS_Montage(name, ra, dec, pix_width, map_width, band, input_dir, out_dir):
#    print('Montaging '+name+'_SDSS_'+band)
#    location_string = str(ra)+' '+str(dec)
#    if os.path.exists(input_dir+'Montage_Temp'):
#        shutil.rmtree(input_dir+'Montage_Temp')
#    os.makedirs(input_dir+'Montage_Temp')
#    os.chdir(input_dir+'Montage_Temp')
#    montage_wrapper.commands.mHdr(location_string, map_width, input_dir+'Montage_Temp/'+str(name)+'_HDR', pix_size=pix_width)
#    montage_wrapper.commands.mExec('SDSS', band.lower(), raw_dir=input_dir, level_only=False, debug_level=0, output_image=out_dir+name+'_SDSS_'+band+'.fits', region_header=input_dir+'Montage_Temp/'+str(name)+'_HDR', workspace_dir=input_dir+'Montage_Temp')
#    shutil.rmtree(input_dir+'Montage_Temp')
#    print('Completed Montaging '+name+'_SDSS_'+band)



# Define function to wget SDSS fields; reject and re-acquire fields less than 1MB in size
def SDSS_wget(tile_url, path):
    fitsname = tile_url.split('/')[-1:][0]
    print('Acquiring '+fitsname)
    if os.path.exists(path+fitsname):
        os.remove(path+fitsname)
    else:
        success = False
        while success==False:
            try:
                wget.download(tile_url, out=path+fitsname)
                filesize = os.stat(path+fitsname).st_size
                if float(filesize)<1048576.0:
                    raise NameError('File not large enough')
                print('Successful acquisition of '+fitsname)
                success = True
            except:
                print('Failure! Retrying acquistion of '+fitsname)
                time.sleep(0.1)
                success = False

# Define function to extract SDSS bz2 archives
def SDSS_Extract(filepath):
    print('Decompressing file '+str( filepath.split('/')[-1:][0] ))
    os.system('bzip2 -d '+filepath)

# Define function to Montage together contents of folder
def SDSS_Montage(name, ra, dec, pix_width, map_width, band, input_dir, out_dir):
    print('Montaging '+name+'_SDSS_'+band)
    location_string = str(ra)+' '+str(dec)
    if os.path.exists(input_dir+'Montage_Temp'):
        shutil.rmtree(input_dir+'Montage_Temp')
    os.makedirs(input_dir+'Montage_Temp')
    os.chdir(input_dir+'Montage_Temp')
    montage_wrapper.commands.mHdr(location_string, map_width, input_dir+'Montage_Temp/'+str(name)+'_HDR', pix_size=pix_width)
    montage_wrapper.commands.mExec('SDSS', band.lower(), raw_dir=input_dir, level_only=False, debug_level=0, output_image=out_dir+name+'_SDSS_'+band+'.fits', region_header=input_dir+'Montage_Temp/'+str(name)+'_HDR', workspace_dir=input_dir+'Montage_Temp')
    shutil.rmtree(input_dir+'Montage_Temp')
    print('Completed Montaging '+name+'_SDSS_'+band)





# Commence main task
if __name__ == "__main__":

    # State bands of intereset
    bands = ['u','g','r','i','z']

    # Define paths
    in_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SDSS/Temporary_Files/'
    out_dir = '/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SDSS/Mosaics/'

    # Read in source catalogue
    NESS_cat = np.genfromtxt(dropbox+'Work/Tables/NESS/NESS_Sample.csv', delimiter=',', names=True, dtype=None)
    name_list = NESS_cat['name']

    # Read in list of SDSS DR12 primary fields
    dr12_pri = np.genfromtxt('/home/sarumandata2/spx7cjc/NESS/Ancillary_Data/SDSS/SDSS_DR12_Primary_Fields.dat.gz', names=True)
    dr12_pri = [ str(int(pri['RUN']))+' '+str(int(pri['CAMCOL']))+' '+str(int(pri['FIELD'])) for pri in dr12_pri ]

    # Identify sources not yet processed
    already_file = 'SDSS_Already_Processed_List.dat'
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
        print('##### THIS IS JUST GRABBING ORIION, K? #####')
        name = 'Orion'#name_list[i].replace(' ','_')
        ra = ra_list[i]
        dec = dec_list[i]
        time_start = time.time()
        width = 0.25
        print('Processing target '+name)

        # Create field processing dirctories (deleting any prior), and set appropriate Python (ie, SWarp) working directory
        tile_dir = in_dir+str(name)+'/'
        if os.path.exists(tile_dir):
            shutil.rmtree(tile_dir)
        os.makedirs(tile_dir)
        os.makedirs(tile_dir+'Raw')
        os.chdir(tile_dir+'Raw')

        # Use Montage to perform query of SDSS DR9 fields (inluding additional 0.14 degrees to search box to account for field size?)
        coverage = True
        sdss_urls = []
        for band in bands:
            if coverage==False:
                continue
            montage_wrapper.commands.mArchiveList('SDSS', band.lower(), str(ra)+' '+str(dec), width, width, tile_dir+'Query_'+band+'.txt')#width*(2.**0.5), width*(2.**0.5)
            # Check if query returned any results
            if os.stat(tile_dir+'Query_'+band+'.txt').st_size==0:
                coverage = False
                continue
            elif sum(1 for line in open(tile_dir+'Query_'+band+'.txt'))<=3:
                coverage = False
                continue

            # If SDSS coverage present, extract field info from query results, removing duplicate entries
            else:
                montage_urls = np.genfromtxt(tile_dir+'Query_'+band+'.txt', skip_header=2, usecols=[9], dtype=[('URL','S256')])
                [ sdss_urls.append(montage_url[0]) for montage_url in montage_urls ]

        # Remove duplicate URLs from URL list, and only use primary runs
        sdss_urls = list(set(sdss_urls))
        sdss_urls_pri = SDSS_Primary_Check(sdss_urls, dr12_pri)
        print(str(len(sdss_urls_pri))+' of '+str(len(sdss_urls))+' urls for '+name+' are primary')
        sdss_urls = sdss_urls_pri
        if len(sdss_urls)==0:
            coverage = False

        # If no SDSS coverge, record, clean up, and skip onwards
        if coverage==False:
            print(name+' not covered by SDSS')
            alrady_processed_file = open(already_file, 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()
            shutil.rmtree( os.path.join( in_dir, name ) )
            time_list.append(time.time())
            continue
        else:

            # In parallel, download SDSS fields (NB: SDSS server will not permit more than 5 simultaneous wget downloads)
            print('Downloading '+str(len(sdss_urls))+' SDSS fields for '+name)
            dl_complete = False
            dl_fail_count = 0
            while not dl_complete:
                try:
                    signal.alarm(7200)
                    dl_pool = mp.Pool(processes=5)
                    for sdss_url in sdss_urls:
                        input_dir = tile_dir+band+'/'
                        dl_pool.apply_async( SDSS_wget, args=(sdss_url, tile_dir+'Raw/',) )
                    dl_pool.close()
                    dl_pool.join()
                    dl_complete = True
                except:
                    dl_fail_count += 1
                    gc.collect()
                    shutil.rmtree(tile_dir+'Raw')
                    os.makedirs(tile_dir+'Raw')
                    os.chdir(tile_dir+'Raw')
                    if dl_fail_count==5:
                        dl_complete = True
                    print('Download sequence failed; reattemping')

            # In parallel, extract SDSS fields from their bz2 archives
            raw_files = np.array(os.listdir(tile_dir+'Raw/'))
            bz2_pool = mp.Pool(processes=10)
            for raw_file in raw_files:
                bz2_pool.apply_async( SDSS_Extract, args=(tile_dir+'Raw/'+raw_file,) )
            bz2_pool.close()
            bz2_pool.join()

            # Perform metadata query to find if any fields provide coverage at particular coord of target; if not, clean up and move on
            montage_wrapper.commands.mImgtbl(tile_dir+'Raw',  tile_dir+'Metadata_'+band+'.dat', corners=True)
            montage_wrapper.commands_extra.mCoverageCheck(tile_dir+'Metadata_'+band+'.dat', tile_dir+'Metada_Query_'+band+'.dat', ra=ra, dec=dec, mode='point')
            if os.stat(tile_dir+'Metada_Query_'+band+'.dat').st_size>0:
                pass
            elif sum(1 for line in open(tile_dir+'Metada_Query_'+band+'.dat'))>=3:
                pass
            else:
                print(name+' not covered by SDSS')
                no_coverage_file = open('/home/saruman/spx7cjc/NESS/SDSS/SDSS_No_Coverage_List.dat', 'a')
                no_coverage_file.write(name+'\n')
                no_coverage_file.close()
                shutil.rmtree(tile_dir)
                all_complete = True
                continue

            # Copy files to relevant folders
            for band in bands:
                if os.path.exists(tile_dir+band+'/'):
                    shutil.rmtree(tile_dir+band+'/')
                os.makedirs(tile_dir+band+'/')
                raw_files = np.array(os.listdir(tile_dir+'Raw/'))
                for raw_file in raw_files:
                    if '-'+band+'-' in raw_file:
                        shutil.move(tile_dir+'Raw/'+raw_file, tile_dir+band)

            # In parallel, Montage together each band's input fields, whilst dealing with timeouts
            out_pix_width = 0.396
            complete = False
            while not complete:
                try:
                    signal.alarm(3600)
                    pool = mp.Pool(processes=5)
                    for band in bands:
                        input_dir = tile_dir+band+'/'
                        pool.apply_async( SDSS_Montage, args=(name, ra, dec, out_pix_width, width, band, input_dir, out_dir,) )
                        #SDSS_Montage(name, ra, dec, out_pix_width, width, band, input_dir, out_dir)
                    pool.close()
                    pool.join()
                    complete = True
                except:
                    gc.collect()
                    if os.path.exists(os.path.join(input_dir,'Montage_Temp')):
                        shutil.rmtree(os.path.join(input_dir,'Montage_Temp'))
                    print('Mosaicing failed')

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
            print('Estimated completion time: '+time_est)

# Jubilate
print('All done!')


