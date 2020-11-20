# IceMap500
A sea ice detection algorithm based on IceMap250 (Gignac et al., 2017) with new features to generate accurate sea ice extent maps from MODIS visible and infrared data at 500 m resolution at nadir.

***********************************************************************************************************************************

Aiming at the creation of a European sea ice extent indicator, an improved sea ice detection algorithm has been created, mainly based on the classification process used in IceMap250. A new method to correct the effects of NISE artefacts in the MODIS cloud mask and additional threshold tests allow the
enlargement of the mapped area, the reduction of potential error sources and a qualitative improvement of the resulting maps.
Quality assessment has shown this algorithm produces sea ice presence maps systematically achieving accuracies above 90 % (tested in the European Arctic).

This algorithm is distributed as a Python 3 script for Unix systems, requiring commonly used libraries such as numpy and gdal, alongside jenkspy and whitebox libraries (check WhiteboxTools user's manual at https://jblindsay.github.io/ghrg/WhiteboxTools/index.html). Additional software is required to project MODIS swath data and correct the bowtie effect:

-HEG Tool, https://newsroom.gsfc.nasa.gov/sdptoolkit/HEG/HEGHome.html

Inputs are MOD021KM, MOD02HKM and MOD35_L2 hdf files.

***********************************************************************************************************************************

References: Gignac, C., Bernier, M., Chokmani, K., and Poulin, J.: IceMap250- Automatic 250 m Sea Ice Extent Mapping Using MODIS
Data, Remote Sensing, 9, https://doi.org/10.3390/rs9010070, 2017.

Associated datasets: Parera-Portell, J.A., Ubach, R.  IceMap500 European maximum and minimum sea ice extent maps (2000-2019). Universitat Autònoma de
Barcelona,  https://doi.org/10.5565/ddd.uab.cat/233396, 2020.
