# Ancillary-Imaging-Acquisition

A selection of Python scripts for acquiring and standardised imaging data from a wide range of telescopes. Used to assemble ancillary data for the DustPedia, JINGLE, and NESS projects.

They require a number of Python packages to be installed in order to work (including my ChrisFuncs package of convenience functions, available at https://github.com/Stargrazer82301/ChrisFuncs). Many of the scripts need Montage (http://montage.ipac.caltech.edu/), and some use SWarp (https://www.astromatic.net/software/swarp).

Standard workflow to acquire imaginng data for a given telescope is:
 - First, run TelescopeName_Query_Mosaic.py, which will query the online archive for the telescope to find out what data is available, download that data, then mosaic it. 
 - Second, run TelescopeName_Final_Map_Generator.py, which will take the mosaics, calibrate their pixel values to Jy/pix, and generate a standardised header.
 
For both types of script, the user needs to alter script to specify a target catalogue (a csv file with columns giving source name, right ascension, and declination), and input/output directories (plus the path to SWarp, for scripts that use it).

At some point I hope to convert these into a generic, functional form. But at present they are set up to produce ancillary data for the NESS survey, assuming my computer's directory structure.
