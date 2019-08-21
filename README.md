# IceMap500
A sea ice detection algorithm based on IceMap250 (Gignac et al., 2017) with new features to generate accurate sea ice extent maps from MODIS visible and infrared data at 500 m resolution at nadir.

References: Gignac, C., Bernier, M., Chokmani, K., and Poulin, J.: IceMap250- Automatic 250 m Sea Ice Extent Mapping Using MODIS
Data, Remote Sensing, 9, https://doi.org/10.3390/rs9010070, 2017.

Associated datasets: Parera-Portell, J.A., Ubach, R.  European-IceMap250 sea ice extent maps (2000-2017). Universitat Autònoma de
Barcelona,  https://doi.org/10.5565/ddd.uab.cat/196007, 2018.

***********************************************************************************************************************************

Aiming at the creation of a European sea ice extent indicator, an improved sea ice detection algorithm has been created, mainly based on the classification process used in IceMap250. A new method to correct the effects NISE artefacts in the MODIS cloud mask and additional threshold tests allow the
enlargement of the mapped area, the reduction of potential error sources and a qualitative improvement of the resulting maps.
Quality assessment has shown this algorithm produces sea ice presence maps systematically achieving accuracies above 90 % (tested in the European Arctic).

This algorithm is distributed as a Python 3 script adapted for Unix systems, requiring commonly used libraries such as numpy and gdal, alongside the jenkspy library and the WhiteboxTools Python geoprocessing utilities. Additional software is required to project MODIS swath data and correct the bowtie effect:

-HEG Tool, https://newsroom.gsfc.nasa.gov/sdptoolkit/HEG/HEGHome.html

Instructions on how to run the script and required data are detailed in the py file.
Inputs are MOD021KM, MOD02HKM and MOD35_L2 hdf files.
