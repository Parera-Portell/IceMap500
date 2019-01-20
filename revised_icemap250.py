# revised_icemap250.py (24/08/2018)
# Created in: May 2018
# Authors: Joan A. Parera Portell, Raquel Ubach (UAB)
# Original IceMap250 algorithm: Charles Gignac & Jimmy Poulin (INRS, Canada)
# Description: European IceMap250 algorithm
# -----------------------------------------------------------------------------

"""
*****************************ICEMAP250-EUROPE**********************************

This script will generate 250 m resolution sea ice extent maps from MODIS imagery
for the selected marine regions, months and years. 
The core of the script is based on the IceMap250 algorithm (Gignac et al., 2017).
Tested in ArcGIS 10.5 and Python 3.

            USAGE - input variables (must be set below this docstring)

month_list = []          #list of months to be processed as integers
year_list = []           #list of years to be processed as integers
seareg = ""              #water mask (used for land extraction)
a_data = ""              #location of MODIS hdf files: MOD021KM.hdf, 
                          MOD02HKM.hdf, MOD02QKM.hdf, MOD03.hdf, MOD35_L2.hdf
b_bands = ""             #output folder for tif MODIS bands
c_masks = ""             #output folder for MOD35 and VIS masks
d_ndsii2 = ""            #output folder for NDSII-2 files
e_maskeddata = ""        #output folder for masked datasets
f_composite = ""         #output folder for final composite maps
workfolder = ""          #location where temporary files will be stored
ccrs_folder = ""         #location of the CCRS' downscaling executables.
                          Available at ftp://ftp.ccrs.nrcan.gc.ca/ad/
                          CCRS_CANADA/Software/Compositing_V5.5
mrt_folder = ""          #location of NASA's MRTSwath bin directory. Available
                          at https://lpdaac.usgs.gov/tools/modis_reprojection_
                          tool_swath

Projection parameters:
-The output projection is WGS84 Lambert Conformal Conic (spheroid = 6378137.0,
 298.257223563). It cannot be changed as other projections might cause the
 CCRS tool to crash. Bear in mind that to quantify sea ice extent, maps must 
 be in an equal-area projection.
 The following variables must be set by the user (default decimal separator is 
 ".", be sure to change the rest of the projection parameters in the script if 
 your separator in ArcGis is ",")-
    
origin_lat = ""          #latitude of origin as string
origin_lon = ""          #longitude of origin as string
std_parallel_1 = ""      #1st standard parallel of LCC projection as string
std_parallel_2 = ""      #2nd standard parallel of LCC projection as string

*******************************************************************************
"""


#**************************USER-DEFINED VARIABLES******************************
month_list = []
year_list = []
seareg = ""
a_data = ""
b_bands = ""
c_masks = ""
d_ndsii2 = ""
e_maskeddata = ""
f_composite = ""
workfolder = ""
ccrs_folder = ""
mrt_folder = ""

# Projection parameters
origin_lat = ""
origin_lon = ""
std_parallel_1 = ""
std_parallel_2 = ""

#******************************************************************************


# Import modules
import arcpy, os, subprocess, gc, io, numpy
from arcpy.sa import *

# Arcpy environmental variables
arcpy.env.addOutputsToMap = False
arcpy.env.pyramid = "NONE"
arcpy.env.compression = "LZ77"

# Creation of output directories structure
folders = [b_bands, c_masks, d_ndsii2, e_maskeddata, f_composite]
for folder in folders:
    if not os.path.exists(folder):
        os.mkdirs(folder)
        years = year_list
        for yr in years: 
            year = str(folder) + "\\" + str(yr)
            if not os.path.exists(year):
                os.mkdir(year)
                for month in month_list:
                    os.mkdir(str(year) + "\\" + str(month))

# Creation of workfolder
if not os.path.exists(workfolder):
        os.mkdirs(workfolder)


#************************STEP 1: EXTRACTION OF MOD35_L2************************    
# Function to extract binary masks from MOD35_L2
def mask_extraction(mod35, mod35_array, result, bit_position, bit_length, value):
    bit_value = int(value, 2)
    length_str = ""
    for i in range(bit_length):
        length_str += "1"
    bit_length = int(length_str, 2)
    
    value1 = bit_length << bit_position
    value2 = bit_value << bit_position
    array_mask = (mod35_array & value1) == value2
    
    array_int = array_mask.astype(int)
    raster_mask = arcpy.NumPyArrayToRaster(array_int, arcpy.Point(mod35.extent.XMin, mod35.extent.YMin),
                                           mod35, mod35, mod35.noDataValue)
    raster_mask.save(result)


def step1(y,m,d):
    arcpy.AddMessage("STEP 1: EXTRACTION OF MOD35_L2")
    infolder = a_data
    outfolder = c_masks +"\\"+ str(y) +"\\" + m
    arcpy.env.workspace = infolder
    arcpy.env.cellSize = 1000
    mod03 = arcpy.ListDatasets("MOD03.A" + str(y) + str(d) + "*", "Raster")
    if mod03:
        for i in mod03:
            if not os.path.exists(outfolder + "\\" + str(i[6:19]) + "_Cloud_Mask.tif"):
                arcpy.AddMessage(str(len(mod03)) + " MOD03 rasters found.")
                mod35_l2_mask = outfolder + "\\" + str(i[6:19]) + "_Cloud_Mask_b0.tif"
                mod35_l2_m = outfolder + "\\" + str(i[6:19]) + "_Cloud_Mask.tif"
                cloud_binmask = workfolder+"\\"+str(i[6:19]) + "_binmask.tif"
                night_binmask= workfolder+"\\"+str(i[6:19]) + "_nightmask.tif"
                sun_glint = workfolder+"\\"+str(i[6:19]) + "_glint.tif"
                binmask = workfolder+"\\"+str(i[6:19])+ "_binary.tif"
                mod35_l2 = arcpy.ListDatasets("MOD35_L2.A"+str(y)+str(d)+i[14:19]+"*", "Raster")
                
#               Conversion & projection of hdf files with MRTSwath
                if mod35_l2:
                    arcpy.AddMessage("Converting "+ str(i[6:19]))                    
                    input_file = "-if=" + infolder + "\\" + str(mod35_l2[0])
                    input_geolocation= "-gf=" + infolder + "\\" + str(i)
                    output_file = "-of=" + outfolder + "\\" + str(i[6:19])
                    dataset = "-sds=Cloud_Mask,1,0,0,0,0,0"
                    output_format = "-off=GEOTIFF_FMT"
                    kernel = "-kk=NN"
                    projection= "-oproj=PS"
                    projection_param = "-oprm=0.0,0.0,0.0,0.0,0.0,90.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0"
                    pixel_size = "-opsz=1000.0"
                    sphere = "-osp=8"
                    path_mrt = mrt_folder + "\\swath2grid"         
                    subprocess.call([path_mrt, input_file, input_geolocation, output_file, dataset, output_format, kernel, projection, projection_param, pixel_size, sphere])
                
                    if not os.path.exists(binmask) and not os.path.exists(mod35_l2_m):
                        arcpy.DefineProjection_management(mod35_l2_mask, "PROJCS['North_Pole_Stereographic',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Stereographic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',90.0],UNIT['Meter',1.0]]")
                        arcpy.ProjectRaster_management(mod35_l2_mask, mod35_l2_m, "PROJCS['LC_WGS84',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',"+origin_lon+"],PARAMETER['Standard_Parallel_1',"+std_parallel_1+"],PARAMETER['Standard_Parallel_2',"+std_parallel_2+"],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',"+origin_lat+"],UNIT['Meter',1.0]]", "NEAREST", "1000 1000", "", "", "PROJCS['North_Pole_Stereographic',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Stereographic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',90.0],UNIT['Meter',1.0]]")               
                        arcpy.Delete_management(mod35_l2_mask)
                        arcpy.AddMessage("1/3 done")
                              
                        # Extraction of cloud mask Quality Assessment Bits
                        mod35_array = arcpy.RasterToNumPyArray(mod35_l2_m)
                        if not os.path.exists(cloud_binmask):
                            mask_extraction(mod35_l2_m, mod35_array, cloud_binmask, 1, 2, "11")
                            arcpy.DefineProjection_management(cloud_binmask, "PROJCS['LC_WGS84',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',"+origin_lon+"],PARAMETER['Standard_Parallel_1',"+std_parallel_1+"],PARAMETER['Standard_Parallel_2',"+std_parallel_2+"],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',"+origin_lat+"],UNIT['Meter',1.0]]")
                        if not os.path.exists(night_binmask):
                            mask_extraction(mod35_l2_m, mod35_array, night_binmask, 3, 1, "1")
                            arcpy.DefineProjection_management(night_binmask, "PROJCS['LC_WGS84',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',"+origin_lon+"],PARAMETER['Standard_Parallel_1',"+std_parallel_1+"],PARAMETER['Standard_Parallel_2',"+std_parallel_2+"],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',"+origin_lat+"],UNIT['Meter',1.0]]")
                        if not os.path.exists(sun_glint):
                            mask_extraction(mod35_l2_m, mod35_array, sun_glint, 4, 1, "1")
                            arcpy.DefineProjection_management(sun_glint, "PROJCS['LC_WGS84',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',"+origin_lon+"],PARAMETER['Standard_Parallel_1',"+std_parallel_1+"],PARAMETER['Standard_Parallel_2',"+std_parallel_2+"],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',"+origin_lat+"],UNIT['Meter',1.0]]")
                        arcpy.AddMessage("2/3 done")
                        arcpy.Delete_management(mod35_array)
                              
                        # Setting NoData areas as null
                        if not os.path.exists(binmask):
                            result1 = Con(Raster(mod35_l2_m) != 0, Con(Raster(night_binmask) == 1, Con(Raster(sun_glint) == 1, Raster(cloud_binmask))))
                            result1.save(binmask)    
        gc.collect()
        arcpy.AddMessage("Step 1 completed")
 
       
#**********************STEP 2: CREATION OF MOD35 CLOUD MASK********************
def step2(y,m,d):
    arcpy.AddMessage("STEP 2: CREATION OF MOD35 CLOUD MASK")
    outfolder = c_masks +"\\"+ str(y) +"\\" + m
    arcpy.env.workspace = workfolder
    arcpy.env.snapRaster = seareg
    arcpy.env.cellSize = 250
    binarymasklist = arcpy.ListDatasets("A"+str(y)+str(d)+"*_binary.tif", "Raster")
    if binarymasklist:
        for binmask in binarymasklist:
            mod_sea = workfolder+"\\"+str(binmask[:13]) + "_MOD35seareg.tif"
            cloud_binmask = workfolder+"\\"+str(binmask[:13])+"_binmask.tif"
            night_binmask = workfolder+"\\"+str(binmask[:13])+"_nightmask.tif"
            sun_glint = workfolder+"\\"+str(binmask[:13])+"_glint.tif"
            mod35 = outfolder +"\\"+ str(binmask[:13]) + "_MOD35.tif"
            
            # Extraction of cloud mask by Sea Regions
            if not os.path.exists(mod35):
                try:
                    arcpy.AddMessage("Processing"+str(binmask[:13]))
                    result1 = Con(Raster(seareg) == 1, Raster(binmask))
                    result1.save(mod_sea)
                    if arcpy.GetRasterProperties_management(mod_sea, "ALLNODATA").getOutput(0) == "1":
                        arcpy.Delete_management([cloud_binmask, night_binmask, sun_glint, binmask, mod_sea])
                        arcpy.AddMessage("MOD35 cloud mask out of range!")
                    else:
                        result2 = 1 * Raster(mod_sea)
                        result2.save(mod35)
                        arcpy.AddMessage("MOD35 cloud mask created!") 
                        arcpy.Delete_management([cloud_binmask, night_binmask, sun_glint, binmask, mod_sea])
                except:
                    arcpy.Delete_management([cloud_binmask, night_binmask, sun_glint, binmask])
                    arcpy.AddMessage("MOD35 cloud mask out of range!")
    gc.collect()
    arcpy.AddMessage("Step 2 completed")
 
      
#**********************STEP 3: DOWNSCALING OF BANDS 4 & 7**********************
def step3(y,m,d):
    arcpy.AddMessage("STEP 3: DOWNSCALING OF BANDS 4 & 7")
    infolder = a_data
    arcpy.env.workspace = infolder
    mod02qkm = arcpy.ListDatasets("MOD02QKM.A"+str(y)+str(d)+"*", "Raster")
    if mod02qkm:
        for modqkm in mod02qkm:
            if not os.path.exists(infolder+"\\MOD02DKM.A"+modqkm[10:]):
                if os.path.exists(c_masks+"\\"+str(y)+"\\"+ m +"\\"+str(modqkm[9:22])+"_MOD35.tif"):
                    arcpy.AddMessage(str(len(mod02qkm)) + " MOD02QKM rasters found.")
                    mod02hkm = arcpy.ListDatasets("MOD02HKM.A"+str(y)+str(d)+modqkm[17:22]+"*", "Raster")
                    if mod02hkm:
                        arcpy.AddMessage("Processing" + str(modqkm[9:]))
                        f = open(infolder + "\\downsc_List.txt", 'w')
                        f.write(str(infolder)+ "\\" + str(mod02hkm[0]) + "\n"
                                + str(infolder)+ "\\" + str(modqkm))
                        f.close()
                        path_crms = ccrs_folder +"\\"+ 'mod02_500m_to_250m_DownScale_v3.0.exe'
                        subprocess.call([path_crms, "path="+infolder, "file=downsc_List.txt", "band=4,7"])
                        arcpy.AddMessage(str(mod02hkm[0])+" downscaled!")
        arcpy.AddMessage("Step 3 completed")


#********STEP 4: REPROJECTION AND CONVERSION OF BANDS 02, 04, 20 & 32**********
def step4(y,m,d):
    arcpy.AddMessage("STEP 4: REPROJECTION AND CONVERSION OF BANDS 2, 4, 20 & 32")
    infolder = a_data
    arcpy.env.workspace = infolder
    mod02dkm = arcpy.ListDatasets("MOD02DKM.A"+str(y)+str(d)+"*", "Raster")
    if mod02dkm:
        for moddkm in mod02dkm:
            if not os.path.exists(b_bands + "\\"+ str(y) +"\\"+ m + "\\" + str(moddkm[9:22]) + "_B32_final.tif"):
                if os.path.exists(c_masks+"\\"+str(y)+"\\"+ m +"\\"+str(moddkm[9:22])+"_MOD35.tif"):
                    mod02qkm = arcpy.ListDatasets("MOD02QKM.A"+str(y)+str(d)+moddkm[17:22]+"*", "Raster")
                    mod021km = arcpy.ListDatasets("MOD021KM.A"+str(y)+str(d)+moddkm[17:22]+"*", "Raster")
                    mod03 = arcpy.ListDatasets("MOD03.A"+str(y)+str(d)+moddkm[17:22]+"*", "Raster")
                    if mod02qkm and mod021km and mod03:
                        arcpy.AddMessage("Processing" + str(mod03[0][6:19]))                    
                        # Creating configuration file for 250m resolution data
                        f = io.open(infolder + '\\CMRP_250.cfg', 'w')
                        f.write("Version\n" + "2.5\n"
                                + "CfgFile_SeparateSpaces\n" + "1\n" 
                                + "Projection_Type\n"  + "LCC\n"
                                + "E0\n" + "4321000\n" 
                                + "N0\n"+ "3210000\n" 
                                + "E0_std\n" + "4321000\n"
                                + "N0_std\n" + "3210000\n"
                                + "pix_size\n" + "250\n"
                                + "Lat_Origin\n" + origin_lat + "\n"
                                + "Lon_Origin\n" + origin_lon + "\n"
                                + "LCC_Lat1\n" + std_parallel_1 + "\n"
                                + "LCC_Lat2\n" + std_parallel_1 + "\n"
                                + "Width\n" + "0\n"
                                + "Height\n" + "0\n"
                                + "Auto_Reduce\n" + "1\n"
                                + "Save_Granule_Shape\n" + "0\n"
                                + "Smooth_LatLon\n"  + "0\n"
                                + "Offset_Boundary\n" + "0\n"
                                + "1km_Channels\n" + "-\n"
                                + "500m_Channels\n" + "-\n"
                                + "250m_Channels\n" + "2 4 7\n"
                                + "Reproject_Angles\n" + "0\n"
                                + "MOD03_File_Name\n" + str(mod03[0]) + "\n"
                                # MOD021KM & MOD02HKM are not inputs
                                + "MOD021_File_Name\n"  + "-\n"
                                + "MOD02H_File_Name\n"  + "-\n"
                                # end of unused inputs
                                + "MOD02Q_File_Name\n" + str(mod02qkm[0]) + "\n"
                                + "MOD02DS_File_Name\n"  + str(moddkm) + "\n"
                                + "Output_File_Name_Template\n"
                                + workfolder + "\\" + str(mod03[0][6:19]) + "_B%n" + ".raw\n"
                                + "Reproject_BiLinear\n" + "1\n"
                                + "Convert_to_Reflectance\n" + "1\n"
                                + "Convert_to_Temperature\n" + "0\n"
                                + "Correct_Low_Saturation\n" + "0\n"
                                + "Correct_High_Saturation\n" + "0\n"
                                + "Apply_Filtering\n" + "0\n"
                                + "Write_AuxFile\n" + "1\n"
                                + "ScaleOut_Ref\n" + "1000\n"
                                + "ScaleOut_Rad\n" + "100.\n"
                                + "ScaleOut_Temp\n" + "100")
                        f.close()
                        
                        # Creating configuration file for 1km resolution data
                        f = open(infolder+'\\CMRP_1000.cfg', 'w')
                        f.write("Version\n" + "2.5\n"
                                + "CfgFile_SeparateSpaces\n" + "1\n"
                                + "Projection_Type\n" + "LCC\n"
                                + "E0\n"  + "4321000\n" 
                                + "N0\n" + "3210000\n" 
                                + "E0_std\n" + "4321000\n"
                                + "N0_std\n" + "3210000\n"
                                + "pix_size\n" + "1000\n"
                                + "Lat_Origin\n" + origin_lat + "\n"
                                + "Lon_Origin\n" + origin_lon + "\n"
                                + "LCC_Lat1\n" + std_parallel_1 + "\n"
                                + "LCC_Lat2\n" + std_parallel_1 + "\n"
                                + "Width\n" + "0\n"
                                + "Height\n" + "0\n"
                                + "Auto_Reduce\n" + "1\n"
                                + "Save_Granule_Shape\n" + "0\n"
                                + "Smooth_LatLon\n" + "0\n"
                                + "Offset_Boundary\n" + "0\n"
                                + "1km_Channels\n" + "20 32\n"
                                + "500m_Channels\n" + "-\n"
                                + "250m_Channels\n" + "-\n"
                                + "Reproject_Angles\n" + "0\n"
                                + "MOD03_File_Name\n" + str(mod03[0]) + "\n"
                                + "MOD021_File_Name\n" + str(mod021km[0]) + "\n"
                                # MOD02HKM, MOD02QKM & MOD02DKM are not inputs
                                + "MOD02H_File_Name\n" + "-\n"
                                + "MOD02Q_File_Name\n" + "-\n"
                                + "MOD02DS_File_Name\n" + "-\n"
                                # end of unused inputs
                                + "Output_File_Name_Template\n"
                                + workfolder + "\\" + str(mod03[0][6:19]) + "_B%n" + ".raw\n"
                                + "Reproject_BiLinear\n" + "1\n"
                                + "Convert_to_Reflectance\n" + "0\n"
                                + "Convert_to_Temperature\n" + "1\n"
                                + "Correct_Low_Saturation\n" + "0\n"
                                + "Correct_High_Saturation\n" + "0\n"
                                + "Apply_Filtering\n" + "0\n"
                                + "Write_AuxFile\n" + "1\n"
                                + "ScaleOut_Ref\n" + "1000.\n"
                                + "ScaleOut_Rad\n" + "1000\n"
                                + "ScaleOut_Temp\n" + "100.")
                        f.close()
                        
                        path_cmrp = ccrs_folder + "\\" + 'CMRP_v257.exe'
                        subprocess.call([path_cmrp, "-v2", "-c"+infolder+"\\CMRP_250.cfg", infolder])
                        subprocess.call([path_cmrp, "-v2", "-c"+infolder+"\\CMRP_1000.cfg", infolder])
                        arcpy.AddMessage(str(mod03[0][6:19])+" data processed!")
    arcpy.AddMessage("Step 4 completed")
   
    
#**************************STEP 5: FORMAT CONVERSIONS**************************
def step5(y,m,d):
    arcpy.AddMessage("STEP 5: FORMAT CONVERSIONS")
    outfolder = b_bands +"\\"+ str(y) +"\\"+ m
    maskfolder = c_masks +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = workfolder
    arcpy.env.snapRaster = seareg
    arcpy.env.cellSize = 250
    bandlist = arcpy.ListDatasets("A"+str(y)+str(d)+"*.raw", "Raster")
    if bandlist:
        for band in bandlist:
            band_1 = workfolder + "\\" + str(band[:17]) + ".tif"
            band_clip = workfolder + "\\" + str(band[:17]) + "_clip.tif"
            band_final= outfolder +"\\"+ str(band[:17]) + "_final.tif"
            
            if not os.path.exists(band_final) and os.path.exists(maskfolder +"\\" + str(band[:13]) + "_MOD35.tif"):
                arcpy.AddMessage("Processing "+str(band)+"...")
                
                # Convertion bands from raw to TIFF
                if not os.path.exists(band_1):
                    arcpy.CopyRaster_management(band, band_1, "", "", "65535", "NONE", "NONE", "16_BIT_UNSIGNED", "NONE", "NONE", "TIFF", "NONE")
                arcpy.AddMessage("1/5 done")
                
                # Projection definition
                arcpy.DefineProjection_management(band_1, "PROJCS['LC_WGS84',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',"+origin_lon+"],PARAMETER['Standard_Parallel_1',"+std_parallel_1+"],PARAMETER['Standard_Parallel_2',"+std_parallel_2+"],PARAMETER['Scale_Factor',1.0],PARAMETER['Latitude_Of_Origin',"+origin_lat+"],UNIT['Meter',1.0]]")
                arcpy.AddMessage("2/5 done")
                
                # Extraction by Sea Regions mask
                if not os.path.exists(band_clip):
                    try:
                        result = Con(Raster(seareg) == 1, Raster(band_1))
                        result.save(band_clip)
                        arcpy.AddMessage("3/5 done")
                        
                        # Setting NoData areas as null
                        if not os.path.exists(band_final):
                            if arcpy.GetRasterProperties_management(band_clip, "MAXIMUM").getOutput(0) == "0":
                                arcpy.Delete_management([band, band_1, band_clip])
                                arcpy.AddMessage("Scene out of range")
                                gc.collect()
                            else:
                                arcpy.gp.SetNull_sa(band_clip, band_clip, band_final, "VALUE = 0")
                                arcpy.Delete_management([band, band_1, band_clip])
                                arcpy.AddMessage("Band converted to tif!")
                                gc.collect()
                    except:
                        arcpy.Delete_management([band, band_1])
                        arcpy.AddMessage("Scene out of range")
    gc.collect()
    arcpy.RefreshCatalog(workfolder)
    arcpy.RefreshCatalog(outfolder)
    arcpy.AddMessage("Step 5 completed")


#************************STEP 6: GENERATION OF VIS MASK************************
def step6(y,m,d):
    arcpy.AddMessage("STEP 6: GENERATION OF VIS MASK")
    arcpy.AddMessage("Generating VIS mask...")
    infolder = b_bands +"\\"+ str(y) +"\\"+ m
    outfolder = c_masks +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = infolder
    band20list = arcpy.ListDatasets("A"+str(y)+str(d)+"*20_final*", "Raster")
    if band20list:
        for band20 in band20list:
            if not os.path.exists(outfolder +"\\"+ str(band20[:13]) + "_VISmask.tif"):
                arcpy.AddMessage("Processing "+str(band20[:13])+"...")
                band32 = str(band20[:15]) + "32_final.tif"
                r20_32 = workfolder + "\\" + str(band20[:13]) + "_R.tif"
                vis = workfolder + "\\" + str(band20[:13]) + "_VISeq.tif"
                mask_vis = outfolder +"\\"+ str(band20[:13]) + "_VISmask.tif"
                
                # Calculation of normalized difference between bands 20 and 32
                if not os.path.exists(r20_32):
                    result_0 = Int(((Raster(band20)-Raster(band32))*10000)/(Raster(band20)+Raster(band32)))
                    result_0.save(r20_32)
                arcpy.AddMessage("1/4 done")
        
                # Get Raster Statistics: mean and std
                mean_val = arcpy.GetRasterProperties_management(r20_32, "MEAN").getOutput(0)
                mean_val = mean_val.replace(",",".")
                mean_fl = float(mean_val)
                std = arcpy.GetRasterProperties_management(r20_32, "STD").getOutput(0)
                std = std.replace(",",".")
                std_fl = float(std)
                arcpy.AddMessage("2/4 done")
        
                # Calculation of VIS equation
                if not os.path.exists(vis):
                    result_1 = Int(((Raster(r20_32)-(mean_fl))*1000)/ std_fl)
                    result_1.save(vis)
                arcpy.AddMessage("3/4 done")
        
                # Threshold VIS mask
                if not os.path.exists(mask_vis):
                    result_2 = Con(Raster(vis)<500,1,0)
                    result_2.save(mask_vis)
                arcpy.Delete_management([r20_32, vis])
                arcpy.AddMessage("VIS mask created!") 
        gc.collect()
        arcpy.RefreshCatalog(workfolder)
        arcpy.RefreshCatalog(outfolder)
        arcpy.AddMessage("Step 6 completed")  


#*************************STEP 7: CALCULATION OF NDSII-2***********************
def step7(y,m,d):
    arcpy.AddMessage("STEP 7: CALCULATION OF NDSII-2")
    arcpy.AddMessage("Calculating NDSII-2...")
    infolder = b_bands +"\\"+ str(y) +"\\"+ m
    outfolder = d_ndsii2 +"\\"+ str(y) +"\\"+ m
    modfolder = c_masks +"\\"+ str(y) +"\\"+ m
    visfolder = c_masks +"\\"+ str(y) +"\\"+ m
    maskdatafolder = e_maskeddata +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = infolder
    band02list = arcpy.ListDatasets("A"+str(y)+str(d)+"*02_final*", "Raster")
    if band02list:
        for band02 in band02list:
            arcpy.AddMessage("Processing "+str(band02[:13])+"...")
            band04 = infolder + "\\" + str(band02[:15]) + "04_final.tif"
            mask_mod = modfolder +"\\"+ str(band02[:13]) + "_MOD35.tif"
            mask_vis = visfolder +"\\"+ str(band02[:13]) + "_VISmask.tif"
            green_modm = workfolder + "\\" + str(band02[:15]) + "04_m.tif"
            swir_modm = workfolder + "\\" + str(band02[:17]) + "_m.tif"
            ndsii2 = outfolder +"\\"+ str(band02[:13]) + "_NDSII2.tif"
            
            if not os.path.exists(maskdatafolder +"\\"+str(band02[:13]) + "_MOD35map.tif") and not os.path.exists(ndsii2):
                if os.path.exists(mask_mod) and os.path.exists(mask_vis) and os.path.exists(band04):
                # Deletion of out of range values
                    if not os.path.exists(green_modm):
                        result_0 = Con(Raster(band04) < 1001, Raster(band04))
                        result_0.save(green_modm)
                    if not os.path.exists(swir_modm):
                        result_1 = Con(Raster(band02) < 1001, Raster(band02))
                        result_1.save(swir_modm)
                    arcpy.AddMessage("1/2 completed!")
                    
                # Calculation of NDSII-2
                    if arcpy.GetRasterProperties_management(swir_modm, "ALLNODATA").getOutput(0) == "0" and arcpy.GetRasterProperties_management(green_modm, "ALLNODATA").getOutput(0) == "0":
                        result_2 = Int((((Raster(green_modm) - Raster(swir_modm))*1000) / (Raster(green_modm) + Raster(swir_modm)))+1000)
                        result_2.save(ndsii2)
                        arcpy.Delete_management(swir_modm)
                        arcpy.AddMessage("NDSII-2 calculated!")
                    else:
                        arcpy.Delete_management([swir_modm, green_modm])
                        arcpy.AddMessage("NDSII-2 calculation failed.")
        gc.collect()
        arcpy.RefreshCatalog(workfolder)
        arcpy.RefreshCatalog(outfolder)
        arcpy.AddMessage("Step 7 completed")


#*************************STEP 8: CREATION OF VIS MAP**************************
def step8(y,m,d):
    arcpy.AddMessage("STEP 8: CREATION OF VIS MAP")
    arcpy.AddMessage("Creating VIS map...")
    infolder = d_ndsii2 +"\\"+ str(y) +"\\"+ m
    outfolder = e_maskeddata +"\\"+ str(y) +"\\"+ m
    modfolder = c_masks +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = infolder
    ndsiilist = arcpy.ListDatasets("A"+str(y)+str(d)+"*_NDSII2.tif", "Raster")
    if ndsiilist:
        for index in ndsiilist:
            green = workfolder + "\\" + str(index[:13]) + "_B04_m.tif"
            mask_vis = modfolder +"\\"+ str(index[:13]) + "_VISmask.tif"
            index_vis = workfolder + "\\" + str(index[:13]) + "_indexVIS.tif"
            jenks_vis = workfolder + "\\" + str(index[:13]) + "_jenks.tif"
            map_vis = outfolder +"\\"+ str(index[:13]) + "_VISmap.tif"
            
            if os.path.exists(green):
                # Extract by mask            
                if not os.path.exists(index_vis):
                    arcpy.AddMessage("Processing "+str(index[:13])+"...")
                    result_0 = Con(Raster(mask_vis) == 1,Raster(index))
                    result_0.save(index_vis)
                    if arcpy.GetRasterProperties_management(index_vis, "ALLNODATA").getOutput(0) == "1":
                        arcpy.Delete_management(index_vis)
                    else:
                        # Calculate Jenks natural break
                        if not os.path.exists(jenks_vis):
                            result_1 = Slice(index_vis, "2", "NATURAL_BREAKS", "1")
                            result_1.save(jenks_vis)
                            arcpy.Delete_management(index_vis)
                            arcpy.AddMessage("1/2 done")

                        # Thresholds application
                        if not os.path.exists(map_vis):
                            result_2 = Con(Raster(jenks_vis) == 1, Con(Raster(green) >= 170, 1), Con(Raster(green) < 170, 0, 1))
                            result_2.save(map_vis)     
                arcpy.AddMessage("VIS map created!")
        gc.collect()
        arcpy.RefreshCatalog(workfolder)
        arcpy.RefreshCatalog(outfolder)
        arcpy.AddMessage("Step 8 completed")
  

#*************************STEP 9: CREATION OF MOD35 MAP************************

# STEP 9.0: MOD35 CORRECTION - SCENES WITH ONLY SEA ICE
def step9_0(index, map_mod, b20, b4, slice_b7, outfolder):
    arcpy.AddMessage("Scene with only sea ice. Applying artefact correction.")
    group_mod = workfolder + "\\" + str(index[:13]) + "_MOD35g.tif"
    unique_mod = workfolder + "\\" + str(index[:13]) + "_MOD35u.tif"
    eu_distance = workfolder + "\\" + str(index[:13]) + "_eudist.tif"
    unique_ice = workfolder + "\\" + str(index[:13]) + "_MOD35i.tif"
    mod_c = outfolder + "\\" + str(index[:13]) + "_MOD35c.tif"
    
    result_1 = RegionGroup(Raster(map_mod), "FOUR", "WITHIN", "NO_LINK")
    result_1.save(group_mod)
#   Deletion of sea ice clusters below 1000 pixels
    result_2 = Con(IsNull(Raster(group_mod)), 4, Con(Raster(group_mod), 4, 1, "Count <= 1000"))
    result_2.save(unique_mod)
    if arcpy.GetRasterProperties_management(unique_mod, "MINIMUM").getOutput(0) == "4":
        arcpy.AddMessage("No sea ice in MOD35 correction")
        if not os.path.exists(mod_c):
            result_3 = Con(IsNull(Raster(map_mod)), 4, 4)
            result_3.save(mod_c)
            arcpy.Delete_management([group_mod, unique_mod])
    else:
#       Calculation of euclidean distance from sea ice edge
        if not os.path.exists(eu_distance):
            result_4 = Con(Raster(unique_mod) == 1, 1)
            result_4.save(unique_ice)
            arcpy.Delete_management(unique_mod)
            result_5 = EucDistance(Raster(unique_ice), 35000, 1000)
            result_5.save(eu_distance)
        arcpy.AddMessage("3/4 done")
        if not os.path.exists(mod_c):
#           Sea ice classification in block-affected areas
            arcpy.env.cellSize = 250
            result_6 = Con(((IsNull(Raster(map_mod))) | (IsNull(Raster(unique_ice)))), Con(IsNull(Raster(eu_distance)), 4, Con(IsNull(Raster(slice_b7)), 4, Con(((Raster(b20) < 26600) & (Raster(slice_b7) == 1) & (Raster(b4) >= 170)), 1, 0))), Raster(unique_ice))
            result_6.save(mod_c)
        arcpy.Delete_management([group_mod, unique_ice, eu_distance])
    
    
# STEP 9.1: MOD35 CORRECTION - SCENES WITH ONLY WATER
def step9_1(index, map_mod, outfolder):
    arcpy.AddMessage("Scene with only water. Applying artefact correction.")
    mod_c = outfolder + "\\" + str(index[:13]) + "_MOD35c.tif"
    if not os.path.exists(mod_c):
        result = Con(IsNull(Raster(map_mod)), 4, 0)
        result.save(mod_c)
        
        
# STEP 9.2: MOD35 CORRECTION - SCENES WITH SEA ICE & WATER
def step9_2(index, map_mod, b20, b4, slice_b7, outfolder):
    arcpy.AddMessage("Scene with sea ice & water. Applying artefact correction.")
    group_mod = workfolder + "\\" + str(index[:13]) + "_MOD35g.tif"
    unique_mod = workfolder + "\\" + str(index[:13]) + "_MOD35u.tif"
    eu_distance = workfolder + "\\" + str(index[:13]) + "_eudist.tif"
    unique_ice = workfolder + "\\" + str(index[:13]) + "_MOD35i.tif"
    mod_c = outfolder + "\\" + str(index[:13]) + "_MOD35c.tif"
    
#   Deletion of sea ice clusters below 1000 pixels 
    result_1 = RegionGroup(Raster(map_mod), "FOUR", "WITHIN", "NO_LINK", 0)
    result_1.save(group_mod)  
    result_2 = Con(IsNull(Raster(group_mod)), 4, Con(Raster(group_mod), 4, 1, "Value = 0 OR Count < 1000"))
    result_2.save(unique_mod)
    if arcpy.GetRasterProperties_management(unique_mod, "MINIMUM").getOutput(0) == "4":
        arcpy.AddMessage("No sea ice in MOD35 correction")
        if not os.path.exists(mod_c):
            result_3 = Con(IsNull(Raster(map_mod)), 4, Con(Raster(map_mod) == 1, 4, 0))
            result_3.save(mod_c)
            arcpy.Delete_management([group_mod, unique_mod])
    else:
#       Calculation of euclidean distance from sea ice edge
        if not os.path.exists(eu_distance):
            result_4 = Con(Raster(unique_mod) == 1, 1)
            result_4.save(unique_ice)
            arcpy.Delete_management(unique_mod)
            result_5 = EucDistance(Raster(unique_ice), 35000, 1000)
            result_5.save(eu_distance)
        arcpy.AddMessage("3/4 done")
        if not os.path.exists(mod_c):
#           Sea ice classification in block-affected areas
            arcpy.env.cellSize = 250
            result_6 = Con(((IsNull(Raster(map_mod))) | (IsNull(Raster(unique_ice)))), Con(IsNull(Raster(eu_distance)), 4, Con(IsNull(Raster(slice_b7)), 4, Con(((Raster(b20) < 26600) & (Raster(slice_b7) == 1) & (Raster(b4) >= 170)), 1, 0))), Con(Raster(map_mod) == 0, 0, Con(Raster(unique_ice) == 1, 1, 4)))
            result_6.save(mod_c)
            arcpy.Delete_management([group_mod, unique_ice, eu_distance])
         
     
# STEP 9: GENERATION OF MOD35 MAP and MOD35 ARTEFACT CORRECTION
def step9(y,m,d):
    arcpy.AddMessage("STEP 9: CREATION OF MOD35 MAP")
    arcpy.AddMessage("Creating MOD35 map...")
    infolder = d_ndsii2 +"\\"+ str(y) +"\\"+ m
    outfolder = e_maskeddata +"\\"+ str(y) +"\\"+ m
    modfolder = c_masks +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = infolder
    ndsiilist = arcpy.ListDatasets("A"+str(y)+str(d)+"*_NDSII2.tif", "Raster")
    if ndsiilist:
        for index in ndsiilist:
            green = workfolder +"\\"+ str(index[:13]) + "_B04_m.tif"
            b20 = b_bands +"\\"+ str(y) +"\\"+ m +"\\" + str(index[:13]) + "_B20_final.tif"
            b7 = b_bands +"\\"+ str(y) +"\\"+ m +"\\" + str(index[:13]) + "_B07_final.tif"
            b4 = b_bands +"\\"+ str(y) +"\\"+ m +"\\" + str(index[:13]) + "_B04_final.tif"
            mask_mod = modfolder +"\\"+ str(index[:13]) + "_MOD35.tif"
            map_vis = outfolder +"\\"+ str(index[:13]) + "_VISmap.tif"
            jenks_vis = workfolder +"\\"+ str(index[:13]) + "_jenks.tif"
            map_mod = outfolder +"\\"+ str(index[:13]) + "_MOD35map.tif"
            jenks_mod = workfolder +"\\"+ str(index[:13]) + "_jenksMOD.tif"
            index_mod = workfolder +"\\"+ str(index[:13]) + "_indexMOD.tif"
            index_b7 = workfolder +"\\"+ str(index[:13]) + "_index_b7.tif"
            slice_b7 = workfolder +"\\"+ str(index[:13]) + "_slice_b7.tif"
            
            if not os.path.exists(map_mod):  
#               Check if there is NDSII-2 data after cloud mask application
                arcpy.AddMessage("Processing "+str(index[:13])+"...")
                result_0 = Con(Raster(mask_mod) == 1, Raster(index))
                result_0.save(index_mod)
                if arcpy.GetRasterProperties_management(index_mod, "ALLNODATA").getOutput(0) == "1":
                    arcpy.Delete_management(index_mod)
                    arcpy.AddMessage("NDSII-2 is all NoData.")
                else:
                    # Calculation Jenks natural break 
                    if not os.path.exists(jenks_mod):
                        result_1 = Slice(index_mod, "2", "NATURAL_BREAKS", "1")
                        result_1.save(jenks_mod)
                    if not os.path.exists(index_b7):
                        result_2 = Con(Raster(b7) < 65, Raster(index))
                        result_2.save(index_b7)
                    if arcpy.GetRasterProperties_management(index_b7, "ALLNODATA").getOutput(0) == "0":
                        if not os.path.exists(slice_b7):
                            result_3 = Slice(index_b7, "2", "NATURAL_BREAKS", "1")
                            result_3.save(slice_b7)
                            arcpy.Delete_management([index_b7, index_mod])
                    arcpy.AddMessage("1/4 done")
                    if os.path.exists(jenks_mod) and os.path.exists(map_vis):
                        # MOD35 map threshold tests classification
                        arcpy.env.cellSize = 250
                        result_4 = Con(Raster(jenks_mod) == 1, Con((Raster(green) >= 170) & (Raster(b20) < 26600), 1), Con((Raster(green) < 170), 0))
                        result_4.save(map_mod)
                        arcpy.Delete_management(green)
                        arcpy.AddMessage("2/4 done")
            
            # MOD35 BLOCK CORRECTION
            if os.path.exists(slice_b7) and os.path.exists(map_mod):
                try:
                    # Rasters with only sea ice
                    if arcpy.GetRasterProperties_management(map_mod, "MINIMUM").getOutput(0) == "1":
                        step9_0(index, map_mod, b20, b4, slice_b7, outfolder)                        
                      
                    # Rasters with only water
                    elif arcpy.GetRasterProperties_management(map_mod, "MAXIMUM").getOutput(0) == "0":
                        step9_1(index, map_mod, outfolder)
                    
                    # Rasters with sea ice + water
                    else:
                        step9_2(index, map_mod, b20, b4, slice_b7, outfolder)
                    arcpy.Delete_management([jenks_mod, jenks_vis, slice_b7, index])
                    arcpy.AddMessage("MOD35 map created!")
                    
                except:
                    arcpy.Delete_management([jenks_mod, jenks_vis, slice_b7])
                    arcpy.AddMessage("MOD35 map is all NoData.")
            elif os.path.exists(map_mod):
                mod_c = outfolder + "\\" + str(map_mod[:13]) + "_MOD35c.tif"
                result_5 = Con(IsNull(Raster(map_mod)), 4, Raster(map_mod))
                result_5.save(mod_c)
                arcpy.Delete_management([jenks_mod, jenks_vis, index])
            else:
                arcpy.AddMessage("Scene discarded, no MOD35 map available.")
                arcpy.Delete_management(index)
        gc.collect()
        arcpy.RefreshCatalog(workfolder)
        arcpy.RefreshCatalog(outfolder)
        arcpy.AddMessage("Step 9 completed")
        


#**********************STEP 10: CREATION OF COMPOSITE MAP**********************
def step10(y,m,d):
    arcpy.AddMessage("STEP 10: CREATE COMPOSITE MAP")
    arcpy.AddMessage("Combining maps...")
    infolder = e_maskeddata +"\\"+ str(y) +"\\"+ m
    outfolder = f_composite +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = infolder
    vismaplist = arcpy.ListDatasets("A"+str(y)+str(d)+"*VISmap*", "Raster")
    if vismaplist:
        for vis_m in vismaplist:
            mod_c = infolder + "\\" + str(vis_m[:13]) + "_MOD35c.tif"
            composite_f = workfolder +"\\"+ str(vis_m[:13]) + "_comp.tif"          
            combined_map = workfolder + "\\" + str(vis_m[:13]) + "_combi.tif"
            composite_final = outfolder +"\\"+ str(vis_m[:13]) + "_composite.tif"
            
            if not os.path.exists(composite_final) and not os.path.exists(combined_map):
                arcpy.AddMessage("Processing"+ str(vis_m[:13]))
                # Combination of MOD35 & VIS maps
                if os.path.exists(mod_c):
                    result_0 = (Raster(mod_c) + (Raster(vis_m) *2))
                    result_0.save(combined_map)
                    arcpy.AddMessage("1/2 done")
                    if arcpy.GetRasterProperties_management(combined_map, "MINIMUM").getOutput(0) < "5":
                        # Reclassification of combination output
                        try:
                            remap = RemapValue([[3, 1],[0, 0],[1, 0],[4,0]])
                            result_1 = Reclassify(combined_map, "Value", remap, "NODATA")
                            result_1.save(composite_f)
                            result_2 = Con(Raster(seareg)==1,composite_f)
                            result_2.save(composite_final)
                            arcpy.Delete_management([combined_map, composite_f])
                            arcpy.AddMessage("Compositing done!")
                        except:
                            arcpy.Delete_management(combined_map)
                            arcpy.AddMessage("Compositing failed.")
                    else:
                        arcpy.Delete_management(combined_map)
                        arcpy.AddMessage("Compositing failed. Map is NoData.")
        gc.collect()
        arcpy.RefreshCatalog(workfolder)
        arcpy.RefreshCatalog(outfolder)
        arcpy.AddMessage("Step 10 completed")
   
     
#***********************STEP 11: DAILY MAP SYNTHESIS**********************       
def step11(y,m,d):
    arcpy.AddMessage("STEP 11: SINGLE DAY MAP SYNTHESIS")
    arcpy.AddMessage("Generating map synthesis...")
    infolder = f_composite +"\\"+ str(y) +"\\"+ m
    arcpy.env.workspace = infolder
    day = str(d)+"."
    compositelist = arcpy.ListDatasets("*"+day+"*"+"composite.tif", "Raster")
    if compositelist:
        dailycomposite = "r_" + str(y) + "_" + str(d) + "_daily.tif"
        
        if not os.path.exists(infolder + "\\" + dailycomposite):
                arcpy.AddMessage("Processing day "+ str(d))
                arcpy.env.extent = seareg
                result = CellStatistics(compositelist, "SUM", "DATA")
                result.save(dailycomposite)
    gc.collect()
    arcpy.AddMessage("Step 11 completed")
    
    
#***************STEP 12: MONTHLY MAP SYNTHESIS & EXTENT CALCULATION************
def step12(y,m):
    arcpy.AddMessage("STEP 12: MONTHLY MAP SYNTHESIS")
    arcpy.AddMessage("Generating map synthesis...")
    infolder = f_composite +"\\"+ str(y) +"\\"+ m
    outfolder = f_composite +"\\"+ str(y)
    arcpy.env.workspace = infolder
    arcpy.env.extent = seareg
    complist = arcpy.ListDatasets("*daily*", "Raster")
    synth = infolder + "\\" + m + "_synth.tif"
    synthesis =  outfolder + "\\" + m + "_synthesis.tif"
    seaice_raster = workfolder + "\\" + str(y) +"_"+ m + "_seaice.tif"
    water_raster = workfolder + "\\" + str(y) +"_"+ m + "_water.tif"
    seaice_distance = workfolder + "\\" + str(y) +"_"+ m + "_seaice_dist.tif"
    water_distance = workfolder + "\\" + str(y) +"_"+ m + "_water_dist.tif"
    seaice_ext = workfolder + "\\" + str(y) +"_"+ m + "_extent.tif"
    seaice_extent = outfolder + "\\"+ str(y) +"_"+ m + "_extent.tif"
    
    # Sum of the daily composite maps
    if not os.path.exists(seaice_extent):
        if complist:
            if not os.path.exists(synth):
                result_0 = CellStatistics(complist, "SUM", "DATA")
                result_0.save(synth)
                arcpy.AddMessage("1/5 done")
                
        # Normalization of sea ice map
        if not os.path.exists(synthesis):
            max_val = arcpy.GetRasterProperties_management(synth, "MAXIMUM").getOutput(0)
            max_int = int(max_val)
            result_1 = (Raster(synth)*100)/max_int
            result_1.save(synthesis)
            arcpy.Delete_management(synth)
            arcpy.AddMessage("2/5 done")
        
        # Extraction of sea ice and water
        if not os.path.exists(seaice_raster):
            result_2 = Con(Raster(synthesis) > 5, 1)
            result_2.save(seaice_raster)
            result_3 = Con(Raster(synthesis) == 0, 2)
            result_3.save(water_raster)
            arcpy.AddMessage("3/5 done")
            
        # Euclidean distance calculation
        if not os.path.exists(seaice_distance):
            arcpy.env.cellSize = 250
            result_4 = EucDistance(Raster(seaice_raster))
            result_4.save(seaice_distance)
            result_5 = EucDistance(Raster(water_raster))
            result_5.save(water_distance)
            arcpy.AddMessage("4/5 done")
        
        # Extent calculation
            result_6 = Con((Raster(seaice_distance) == 0), 1, Con(Raster(seaice_distance) < Raster(water_distance), 1, 0))
            result_6.save(seaice_ext)
            arcpy.env.snapRaster = seareg
            result_7 = Con(Raster(seareg) == 1, seaice_ext)
            result_7.save(seaice_extent)
            arcpy.Delete_management([seaice_raster, water_raster, seaice_distance, water_distance, seaice_ext])    
    gc.collect()
    arcpy.RefreshCatalog(outfolder)
    arcpy.AddMessage("Step 12 completed")


#**************************ALGORITHM INITIALIZATION****************************

def run():
    for d in days_range:
        d = "{:03d}".format(d)
        arcpy.AddMessage("Executing day: " + str(d))
        if not os.path.exists(f_composite +"\\"+ str(y) +"\\"+ m + "\\r_" + str(y) + "_" + str(d) + "_daily.tif"):
            step1(y,m,d)
            step2(y,m,d)
            step3(y,m,d)
            step4(y,m,d)
            step5(y,m,d)
            step6(y,m,d)
            step7(y,m,d)
            step8(y,m,d)
            step9(y,m,d)
            step10(y,m,d)
            step11(y,m,d)
    step12(y,m)
          

month_days = {1:range(1,32),2:range(32,61),3:range(60,92),4:range(91,122),
              5:range(121,153),6:range(152,183),7:range(182,214),
              8:range(213,245),9:range(244,275),10:range(274,306),
              11:range(305,336),12:range(335,367)}

for y in year_list:
    for month in month_list:
        days_range = month_days[month]
        m = str(month)
        run()

