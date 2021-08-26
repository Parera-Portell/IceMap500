# IceMap500

A sea ice detection algorithm to map sea ice using MODIS visible and infrared data at 500 m resolution. Tested under Unix systems.

Credits: Joan A. Parera-Portell, Raquel Ubach and Charles Gignac.
See https://doi.org/10.5194/tc-15-2803-2021 for technical details on the algorithm.
***********************************************************************************************************************************

REQUIREMENTS

-Python 3

  Including the following libraries:
  
    -numpy
    
    -osgeo
    
    -jenkspy
    
    -whitebox

-HEG Tool, https://newsroom.gsfc.nasa.gov/sdptoolkit/HEG/HEGHome.html


INPUT FILES

    -MOD021KM / MYD021KM

    -MOD02HKM / MYD02HKM

    -MOD35_L2 / MYD35_L2
  
  
MODIS original hdf data must be stored in a directory structured by year and month:
  
    -->data_dir
  
        -->yyyy
      
            -->mm

All output maps are projected in Lambert Azimuthal Equal Area.

To run the script simply type *python3 [path]/icemap500.py* in your terminal. The script will automatically ask the
user for input and options.
