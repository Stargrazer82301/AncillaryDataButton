# The Ancillary Data Button

A selection of Python scripts for acquiring and standardised imaging data from a wide range of telescopes. Used to assemble ancillary data for the DustPedia (Clark et al., 2018), JINGLE (Saintonge et al., 2018), and NESS projects. If you use these scripts, please cite Clark et al. (2018), which describes their operation in some detail.

They require a number of Python packages to be installed in order to work (including my ChrisFuncs package of convenience functions, available at https://github.com/Stargrazer82301/ChrisFuncs). Many of the scripts need Montage (http://montage.ipac.caltech.edu/), and some use SWarp (https://www.astromatic.net/software/swarp).

I'm currently updating this reposotory, so that a standardised functional form is used for acquiring imaging data for every telescope. Each telescope will have a module called TelescopeName_Button.py. Within that module will be a function called Run, that will require arguments 'ra', 'dec', and 'width'. Calling the function will produce an imagine cutout centred at [ra, dec], with the required width. Alternatively, a sequence of values (list, array, etc) can be provided for these arguments in order to process a list of targets; the function will loop over these targets in turn. Additionally, a number of keyword arguments are available to change how the imaging data is produced. Full details are given in the docstrings of the Run function for each telescope 

So far only some telescopes have had their scripts updated to this new format. The older scripts (found in the NESS folder) follow this workflow:
 - First, run TelescopeName_Query_Mosaic.py, which will query the online archive for the telescope to find out what data is available, download that data, then mosaic it. 
 - Second, run TelescopeName_Final_Map_Generator.py, which will take the mosaics, calibrate their pixel values to Jy/pix, and generate a standardised header. 
 
For both types of script, the user needs to alter script to specify a target catalogue (a csv file with columns giving source name, right ascension, and declination), and input/output directories (plus the path to SWarp, for scripts that use it). The scripts in the NESS folder are set up to produce ancillary data for the NESS survey, assuming my computer's directory structure.
