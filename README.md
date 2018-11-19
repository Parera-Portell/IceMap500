# revised-icemap250
A modified version of IceMap250 (Gignac et al., 2017) with new features to generate accurate sea ice extent maps from MODIS visible and infrared data.

References: Gignac, C., Bernier, M., Chokmani, K., and Poulin, J.: IceMap250- Automatic 250 m Sea Ice Extent Mapping Using MODIS
Data, Remote Sensing, 9, https://doi.org/10.3390/rs9010070, 2017.

Associated datasets: Parera-Portell, J.A., Ubach, R.  European-IceMap250 sea ice extent maps (2000-2017). Universitat Autònoma de
Barcelona,  https://doi.org/10.5565/ddd.uab.cat/196007, 2018.

***********************************************************************************************************************************

Aiming at the creation of a European sea ice extent indicator, a revised version of the previous IceMap250 has been created, using
MODIS downscaled imagery to generate sea ice extent maps at 250 m resolution at nadir. Changes in the classification approach
adapted to the particularities of the European Arctic and a new method to correct NISE artefacts in the MODIS cloud mask allow the
enlargement of the mapped area, the reduction of potential error sources and a qualitative improvement of the resulting maps.
Quality assessment has shown this algorithm produces sea ice presence maps systematically achieving accuracies above 90 %. Such
accuracy levels are not guaranteed when processing other regions.

This algorithm is distributed as a Python script built to work with ArcGIS 10.5. Some additional software is required:

-MODIS Reprojection Tool Swath, https://lpdaac.usgs.gov/tools/modis_reprojection_tool_swath.

-CCRS downscaling and reprojection tools, ftp://ftp.ccrs.nrcan.gc.ca/ad/CCRS_CANADA/Software/Compositing_V5.5.

Instructions on how to run the script and required data are detailed in the py file. It will only work with MODIS L1B data.
