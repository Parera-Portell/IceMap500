# -*- coding: utf-8 -*-
# icemap500.py
# Created in: May 2018
# Author: Joan A. Parera Portell (Universitat Autònoma de Barcelona)

"""
*********************************ICEMAP500*************************************

This script generates improved 500 m resolution sea ice extent maps from MODIS 
imagery.
Based on the IceMap250 algorithm (Gignac et al., 2017).

Maps are projected in Lambert Azimuthal Equal Area, centred at latitudes 90 or 
-90 depending on the hemisphere (EPSGs 102017 and 102020), using the WGS84
ellipsoid.

MODIS hdf data must be in a directory structured by year and month:
    -->hdf_data
        -->2019
            -->03
            
The path to hdf_data will be the directory given to the script as an input. The
intermediate files and results will be stored in a directory following the same
year-month structure.

*******************************************************************************
"""

# Import modules
import numpy as np
import glob, gdal, jenkspy, os, osr, subprocess, sys
from WBT.whitebox_tools import WhiteboxTools
wbt = WhiteboxTools()


#**************************USER-DEFINED VARIABLES******************************
month_list = input("Months (numbers separated by commas): ").split(",")
try:
    month_list = [int(x) for x in month_list]
except:
    print("Months not recognized. Exiting...")
    sys.exit()
    
year_list = input("Years (separated by commas): ").split(",")
try:
    year_list = [int(x) for x in year_list]
except:
    print("Years not recognized. Exiting...")
    sys.exit()
    
hemisphere = input("Hemisphere (N or S): ")
if hemisphere != "N" and hemisphere != "S":
    print("Hemisphere not recognized. Exiting...")
    sys.exit()
    
studyarea = input("Enter path to study area shapefile: ")
if not os.path.exists(studyarea):
    print("Study area shapefile not found. Exiting...")
    sys.exit()
    
heg_bin = input("Enter path to HEG bin folder: ")
if not os.path.exists(heg_bin):
    print("HEG bin not found. Exiting...")
    sys.exit()
 
in_folder = input("Enter input data directory: ")
if not os.path.exists(in_folder):
    print("Input directory not found. Exiting...")
    sys.exit()
    
out_folder = input("Enter directory where data will be stored: ")
if not os.path.exists(out_folder):
    print("Output directory not found. Exiting...")
    sys.exit()


#*********************************FOLDERS*************************************#

hdfs = in_folder+"/data"
bands = out_folder+"/bands"
masks = out_folder+"/masks"
maskedmaps = out_folder+"/maskedmaps"
composites = out_folder+"/composites"
workfolder = out_folder+"/workfolder"

#********************************PROJECTION***********************************#

if hemisphere == "N":
    epsg = 102017
    origin = "90.0"
else:
    epsg = 102020
    origin = "-90.0"

#***************************SUPPLEMENTARY FUNCTIONS***************************#

def directories(year_list, month_list):
    # This function creates the output directory structure.
    folders = [hdfs, bands, masks, maskedmaps, composites]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
        years = year_list
        for yr in years: 
            year = str(folder) + "/" + str(yr)
            if not os.path.exists(year):
                os.mkdir(year)
            for month in month_list:
                monthdir = str(year) + "/" + str("{:02d}".format(month))
                if not os.path.exists(monthdir):
                    os.mkdir(monthdir) 
    # Creation of workfolder
    if not os.path.exists(workfolder):
            os.makedirs(workfolder)

def checkfiles(y,m,d):
    # This function checks whether all files for a given scene are available.
    # Otherwise, files are renamed and excluded from analysis.
    infolder = hdfs+"/"+y+"/"+m
    for mod35hdf in glob.glob(infolder+"/MOD35_L2.A"+y+d+"*"):
        modhkm = glob.glob(infolder+"/MOD02HKM."+mod35hdf[-35:-22]+"*")
        mod1km = glob.glob(infolder+"/MOD021KM."+mod35hdf[-35:-22]+"*")
        if not modhkm or not mod1km:
            scene = mod35hdf.split("/")[-1]
            os.rename(mod35hdf, infolder+"/NA_"+scene)
            if modhkm:
                scene = modhkm[0].split("/")[-1]
                os.rename(modhkm[0], infolder+"/NA_"+scene)
            if mod1km:
                scene = mod1km[0].split("/")[-1]
                os.rename(mod1km[0], infolder+"/NA_"+scene)

def planck(radiance, wvln):
    # This function calculates brightness temperature (wvln in micrometers)
    # using the Inverse Planck Function.
    c1 = 1.191042e8
    c2 = 1.4387752e4
    t = c2/(wvln*np.log(1+(c1/(radiance*wvln**5))))
    return t

def openraster(infile, dtype):
    raster = gdal.Open(infile, gdal.gdalconst.GA_Update)
    array = raster.ReadAsArray().astype(dtype)
    nrows,ncols = np.shape(array)
    trans = raster.GetGeoTransform()
    return array, nrows, ncols, trans, raster

def closeraster(array, outfile, ncols, nrows, dtype, nodata, trans):
    driver = gdal.GetDriverByName("GTiff").Create(outfile, ncols, nrows, 1, 
                                  dtype, options=["COMPRESS=LZW"])
    driver.SetGeoTransform(trans)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    driver.SetProjection(srs.ExportToWkt())
    driver.GetRasterBand(1).SetNoDataValue(nodata)
    driver.GetRasterBand(1).WriteArray(array)
    driver.FlushCache()
    driver = None
    
def hdf_to_tif(infile, outfile, heg_bin, sds, bandnumber, res):
    # This function converts an hdf raster to tif using HEG Tools.
    # Get metadata of hdf file
    raster = gdal.Open(infile)
    band = gdal.Open(raster.GetSubDatasets()[sds][0])
    metadata = band.GetMetadata()
    minx = metadata["WESTBOUNDINGCOORDINATE"]
    maxx = metadata["EASTBOUNDINGCOORDINATE"]
    miny = metadata["SOUTHBOUNDINGCOORDINATE"]
    maxy = metadata["NORTHBOUNDINGCOORDINATE"]
    
    # Create parameter file for HEG
    if res == 500:
        objectname = "MODIS_SWATH_Type_L1B\n"
        fieldname = "EV_500_RefSB|"
    elif res == 250:
        objectname = "MODIS_SWATH_Type_L1B\n"
        fieldname = "EV_250_Aggr500_RefSB|"
    elif res == 1000:
        objectname = "MODIS_SWATH_Type_L1B\n"
        fieldname = "EV_1KM_Emissive|"
    elif res == 0:
        objectname = "mod35"
        fieldname = "Cloud_Mask|"
    else:
        objectname = "mod35"
        fieldname = "Solar_Zenith|"
    
    f = open(heg_bin + "/project.prm", "w")
    f.write("NUM_RUNS = 1\n\n" +
            "BEGIN\n"+
            "INPUT_FILENAME = " + infile + "\n" +
            "OBJECT_NAME = " + objectname + "\n" + 
            "FIELD_NAME = " + fieldname + "\n" +
            "BAND_NUMBER = " + str(bandnumber) + "\n" + 
            "RESAMPLING_TYPE = NN\n" + 
            "OUTPUT_PROJECTION_TYPE = LA\n" +
            "ELLIPSOID_CODE = WGS84\n" +
            "SPATIAL_SUBSET_UL_CORNER = ( " + maxy + " " + minx + " )\n" +
            "SPATIAL_SUBSET_LR_CORNER = ( " + miny + " " + maxx + " )\n" +
            "OUTPUT_PROJECTION_PARAMETERS = ( 0.0 0.0 0.0 0.0 0 "+origin+" 0.0 \
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 )\n" +
            "OUTPUT_PIXEL_SIZE_X = 500.0\n" +
            "OUTPUT_PIXEL_SIZE_Y = 500.0\n" +
            "OUTPUT_FILENAME = " + outfile + "\n" +
            "OUTPUT_TYPE = GEO\n" + 
            "END\n")
    f.close()
    
    # Execute bash commands
    cmd = "LD_LIBRARY_PATH=" + heg_bin + "; export LD_LIBRARY_PATH;\
                    MRTDATADIR=" + heg_bin[:-3] + "data; export MRTDATADIR;\
                    PGSHOME=" + heg_bin[:-3] + "TOOLKIT_MTD; export PGSHOME\
                    ;" + heg_bin + "/swtif -p " + heg_bin + "/project.prm"
    subprocess.run(["/bin/bash", "-c", cmd], check=True)
    
    # Convert SI to reflectance*cos(solar zenith)
    if res == 500 or res == 250:
        ref_scales = metadata["reflectance_scales"].split(",")
        ref_offsets = metadata["reflectance_offsets"].split(",")
        nodata = int(metadata["_FillValue"])
        scale = float(ref_scales[bandnumber-1])
        offset = float(ref_offsets[bandnumber-1])
        array, nrows, ncols, trans, _ = openraster(outfile, np.uint16)
        _ = None
        array_ref = np.where(array == nodata, nodata, scale*(array-offset)*10000)
        array_ref = np.round(array_ref, decimals=0)
        closeraster(array_ref, outfile, ncols, nrows, gdal.GDT_UInt16, nodata, 
                    trans)
        print("MOD02HKM converted to TOA reflectance*cos(solar zenith).")
    
    # Convert SI to brightness temperture
    elif res == 1000:
        rad_scales = metadata["radiance_scales"].split(",")
        rad_offsets = metadata["radiance_offsets"].split(",")
        nodata = int(metadata["_FillValue"])
        scale = float(rad_scales[bandnumber-1])
        offset = float(rad_offsets[bandnumber-1])
        array, nrows, ncols, trans, _ = openraster(outfile, np.single)
        _ = None
        formula = scale*(array-offset)
        array_temp = np.where(array == nodata, nodata, formula)
        if bandnumber == 1:
            t = planck(array_temp, 3.75)
        else:
            t = planck(array_temp, 12.02)
        array_temp = np.where(array_temp > 500, nodata, t*100)
        array_temp = np.round(array_temp, decimals=0)
        closeraster(array_temp, outfile, ncols, nrows, gdal.GDT_UInt16, nodata, 
                    trans)
        print("MOD021KM converted to brightness temperature.")
        
    elif res == 0:
        print("MOD35 converted to tif.")
        
    # Create and fill solar zenith raster
    else:
        nodata = int(metadata["_FillValue"])
        array, nrows, ncols, trans, raster = openraster(outfile, np.int16)
        band = raster.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        gdal.FillNodata(targetBand = band, maskBand = None, maxSearchDist = 5,
                        smoothingIterations = 0)
        array = band.ReadAsArray()
        closeraster(array, outfile, ncols, nrows, gdal.GDT_Int16, nodata, trans)
        print("Solar zenith converted to tif.")  
    
def mask_extraction(infile, outfile, bit_position, bit_length, value):
    # This function extracts masks codified as bits in MOD35_L2.
    array, nrows, ncols, trans, _ = openraster(infile, np.byte)  
    _ = None
    # Selection of bits and creation of mask
    bit_value = int(value, 2)   # conversion to integer of base 2 (binary)
    length_str = ""
    for i in range(bit_length):
        length_str += "1"
    bit_length = int(length_str, 2)
    
    value1 = bit_length << bit_position
    value2 = bit_value << bit_position
    array_mask = (array & value1) == value2
    out_array = array_mask.astype(np.byte)
    closeraster(out_array, outfile, ncols, nrows, gdal.GDT_Byte, 255, trans)
    
def clip(infile, outfile, shapefile, nodata, dtype):
    # This function clips a raster according to an input shapefile.
    gdal.Warp(outfile, infile, format = "GTiff", srcNodata = nodata,
                          dstNodata = nodata, cutlineDSName = shapefile, 
                          cropToCutline = True, outputType = dtype,
                          xRes = 500, yRes = 500, dstSRS = "epsg:"+str(epsg),
                          creationOptions=["COMPRESS=LZW"])
    
def to_ref(infile, outfile, zenith):
    # This function converts reflectance*cos(solar zenith) to reflectance.
    array, nrows, ncols, trans, _ = openraster(infile, np.uint16)
    _ = None
    rows, cols = np.shape(zenith)
    if rows == nrows and cols == ncols:
        con = (zenith != -32767) & (array != 65535)
        rads = np.deg2rad(zenith/100)
        ref = np.where(con, array/np.cos(rads), 65535)
        con2 = con & (ref > 10000)
        ref1 = np.where(con2, 65535, ref)
        ref1 = np.round(ref1, decimals=0)
        closeraster(ref1, outfile, ncols, nrows, gdal.GDT_UInt16, 65535, trans)
        os.remove(infile)
    else:
        print("Rasters are not the same size.")




#**************************STEP1: MOD35 MASK CREATION*************************#

def step1(y,m,d):   
    infolder = hdfs+"/"+y+"/"+m
    checkfiles(y,m,d)
    for mod35hdf in glob.glob(infolder+"/MOD35_L2.A"+y+d+"*"):
        mod35 = bands+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".MOD35.tif"
        index = workfolder+"/index.shp"
        seareg = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".seareg_clip.shp"
        cloud = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".clouds.tif"
        night = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".night.tif"
        glint = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".glint.tif"
        sea = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".sea.tif"
        szenith = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".solar_zenith.tif"
        mask = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".mod35_not_clipped.tif"
        finalmask = masks+"/"+y+"/"+m+"/"+mod35hdf[-35:-22]+".mod35.tif"
        
        # Extract binary masks from the MO35_L2 product
        if not os.path.exists(finalmask):
            if not os.path.exists(mod35):
                hdf_to_tif(mod35hdf, mod35, heg_bin, 7, 1, 0)
            if not os.path.exists(cloud):
                mask_extraction(mod35, cloud, 1, 2, "11")
            if not os.path.exists(night):
                mask_extraction(mod35, night, 3, 1, "1")
            if not os.path.exists(glint):
                mask_extraction(mod35, glint, 4, 1, "1")
            if not os.path.exists(sea):
                mask_extraction(mod35, sea, 6, 2, "00")
            if not os.path.exists(szenith):
                hdf_to_tif(mod35hdf, szenith, heg_bin, 11, 1, 1) 
                
            # Create sea regions mask for the given scene
            if not os.path.exists(seareg):
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(epsg)
                proj = srs.ExportToWkt()
                wbt.layer_footprint(mod35, index)
                f = open(index[:-4]+".prj", "w"); f.write(proj); f.close()
                f = open(index[:-4]+".qpj", "w"); f.write(proj[:-28]+"]")
                f.close()
                wbt.clip(index, studyarea, seareg)
                f = open(seareg[:-4]+".prj", "w"); f.write(proj); f.close()
                f = open(seareg[:-4]+".qpj", "w"); f.write(proj[:-28]+"]")
                f.close()
                  
            # Create MOD35 mask
            if os.path.exists(seareg):
                cloud_array, nrows, ncols, trans, _ = openraster(cloud, np.uint16)
                night_array, _, _, _, _ = openraster(night, np.uint16)
                glint_array, _, _, _, _ = openraster(glint, np.uint16)
                sea_array, _, _, _, _ = openraster(sea, np.uint16)
                mod35_array, _, _, _, _ = openraster(mod35, np.uint16)
                _ = None
                # Conditional operations to obtain mask
                # 0 = clouds, 1 = clear, 2 = sun glint, 255 = land
                con1 = (cloud_array == 1) & (night_array == 1)
                con2 = (glint_array == 1) & (sea_array == 1)
                con3 = (cloud_array == 0) & (sea_array == 1)
                con4 = con3 & (mod35_array != 0) & (night_array == 1)
                con5 = (sea_array == 1) & (mod35_array != 0) & (night_array == 0)
                con = con1 & con2
                clouds = np.where(con4, 0, 255)
                clouds = np.where(con5, 2, clouds)
                mask_array = np.where(con, 1, clouds)
                closeraster(mask_array, mask, ncols, nrows, gdal.GDT_Byte, 255, 
                            trans)
                # Clip mask with study area
                clip(mask, finalmask, seareg, 255, gdal.GDT_Byte)
                mod35_array, _, _, _, _ = openraster(finalmask, np.uint16)
                _ = None
                # Check if map is all nodata, delete if true
                if mod35_array.min() == 255:
                    os.remove(finalmask)
                    print("\n"+mod35+" falls beyond the study area.")
                    for i in [seareg[:-4]+".prj", seareg[:-4]+".qpj", 
                          seareg[:-4]+".dbf", seareg[:-4]+".shp", 
                          seareg[:-4]+".shx", szenith, szenith+".met", mod35,
                          mod35+".met", mod35hdf]:
                        if os.path.exists(i):
                            os.remove(i)
                else:
                    print("MOD35 mask created!")
                
            # Delete files if map is not within the study area
            else:
                print("\n"+mod35+" falls beyond the study area.")
                for i in [seareg[:-4]+".prj", seareg[:-4]+".qpj", 
                          seareg[:-4]+".dbf", seareg[:-4]+".shp", 
                          seareg[:-4]+".shx", szenith, szenith+".met", mod35,
                          mod35+".met", mod35hdf]:
                    if os.path.exists(i):
                        os.remove(i)
            
            # Delete temporal files
            for i in [cloud, night, glint, sea, mask]:
                if os.path.exists(i):
                    os.remove(i)
            index_list = glob.glob(workfolder+"/index*")
            if len(index_list) != 0:
                for i in index_list:
                    os.remove(i)
    checkfiles(y,m,d)

        
#*****************************STEP2: BAND CREATION****************************#
    
def step2(y,m,d):
    infolder = hdfs+"/"+y+"/"+m
    mod02hkm = glob.glob(infolder+"/MOD02HKM.A"+y+d+"*")
    mod021km = glob.glob(infolder+"/MOD021KM.A"+y+d+"*")
    tmp = workfolder+"/temp.tif"
    
    # Solar reflective bands
    for mod in mod02hkm:
        mask = masks+"/"+y+"/"+m+"/"+mod[-35:-22]+".mod35.tif"
        seareg = masks+"/"+y+"/"+m+"/"+mod[-35:-22]+".seareg_clip.shp"
        szenith = masks+"/"+y+"/"+m+"/"+mod[-35:-22]+".solar_zenith.tif"
        zenith = masks+"/"+y+"/"+m+"/"+mod[-35:-22]+".zenith.tif"
        
        if os.path.exists(mask):
            # Crop solar zenith raster
            if not os.path.exists(zenith):
                clip(szenith, zenith, seareg, -32767, gdal.GDT_Int16)
                os.remove(szenith)
                os.remove(szenith+".met")
       
            # Convert to reflectance independent of solar zenith
            band2 = bands+"/"+y+"/"+m+"/"+mod[-35:-22]+".B2.tif"
            band4 = bands+"/"+y+"/"+m+"/"+mod[-35:-22]+".B4.tif"
            band7 = bands+"/"+y+"/"+m+"/"+mod[-35:-22]+".B7.tif"
            zenith_array, _, _, _, _ = openraster(zenith, np.int16)
            _ = None
                
            if not os.path.exists(band2):
                hdf_to_tif(mod, band2, heg_bin, 2, 2, 250)
                clip(band2, tmp, seareg, 65535, gdal.GDT_UInt16)
                to_ref(tmp, band2, zenith_array)
                
            if not os.path.exists(band4):
                hdf_to_tif(mod, band4, heg_bin, 0, 2, 500)
                clip(band4, tmp, seareg, 65535, gdal.GDT_UInt16)
                to_ref(tmp, band4, zenith_array)
                
            if not os.path.exists(band7):
                hdf_to_tif(mod, band7, heg_bin, 0, 5, 500)
                clip(band7, tmp, seareg, 65535, gdal.GDT_UInt16)
                to_ref(tmp, band7, zenith_array)
            
    # Emissive bands
    for mod in mod021km:
        mask = masks+"/"+y+"/"+m+"/"+mod[-35:-22]+".mod35.tif"
        seareg = masks+"/"+y+"/"+m+"/"+mod[-35:-22]+".seareg_clip.shp"
        if os.path.exists(mask):
            band20 = bands+"/"+y+"/"+m+"/"+mod[-35:-22]+".B20.tif"
            band32 = bands+"/"+y+"/"+m+"/"+mod[-35:-22]+".B32.tif"
            if not os.path.exists(band20):
                hdf_to_tif(mod, tmp, heg_bin, 2, 1, 1000)
                clip(tmp, band20, seareg, 65535, gdal.GDT_UInt16)
                os.remove(tmp)
                
            if not os.path.exists(band32):
                hdf_to_tif(mod, tmp, heg_bin, 2, 12, 1000)
                clip(tmp, band32, seareg, 65535, gdal.GDT_UInt16)
                os.remove(tmp)
                print("Bands created!")
                
        
#*******************************STEP3: VIS MASK*******************************#
                
def step3(y,m,d):
    infolder = bands+"/"+y+"/"+m
    bandlist = glob.glob(infolder+"/A"+y+d+".????.B20.tif")
    
    for band20 in bandlist:
        band32 = infolder+"/"+band20[-21:-8]+".B32.tif"
        mask = masks+"/"+y+"/"+m+"/"+band20[-21:-8]+".mod35.tif"
        vismask = masks+"/"+y+"/"+m+"/"+band20[-21:-8]+".vis.tif"
        
        if not os.path.exists(vismask) and os.path.exists(band32):
            b20_array, nrows, ncols, trans, _ = openraster(band20, np.single)
            b32_array, _, _, _, _ = openraster(band32, np.single)
            m_array, _, _, _, _ = openraster(mask, np.uint16)
            _ = None
            # Calculation of normalized difference
            con = (b32_array != 65535) & (b20_array != 65535) & (m_array <= 1)
            con = con & (b20_array + b32_array != 0)
            formula = (b20_array - b32_array)/(b20_array + b32_array)
            index = np.where(con, formula, np.nan)
            # Calculation of standard score
            mean = np.nanmean(index)
            std = np.nanstd(index)
            score = (index-mean)/std
            vis = np.where(con, score*1000, 65535)
            # Classification
            # 0 = cloud, 1 = clear
            con3 = np.where(vis <= 500, 1, 0)
            vis_array = np.where(vis != 65535, con3, 255)
            closeraster(vis_array, vismask, ncols, nrows, gdal.GDT_Byte, 255, 
                        trans)        
            print("VIS mask created!")


#**************************STEP4: MOD35 CLASSIFICATION************************#        

def step4(y,m,d):
    infolder = bands+"/"+y+"/"+m
    bandlist = glob.glob(infolder+"/A"+y+d+".????.B2.tif")
    
    for band2 in bandlist:
        band4 = infolder+"/"+band2[-20:-7]+".B4.tif"
        band20 = infolder+"/"+band2[-20:-7]+".B20.tif"
        mod35 = masks+"/"+y+"/"+m+"/"+band2[-20:-7]+".mod35.tif"
        mod35map = maskedmaps+"/"+y+"/"+m+"/"+band2[-20:-7]+".mod35map.tif"
        
        if not os.path.exists(mod35map):
            b2_array, nrows, ncols, trans, _ = openraster(band2, np.single)
            b4_array, _, _, _, _ = openraster(band4, np.single)
            b20_array, _, _, _, _ = openraster(band20, np.single)
            m_array, _, _, _, _ = openraster(mod35, np.uint16)
            _ = None
            # Calculation of NDSII-2 and Jenks breaks
            con = (b2_array != 65535) & (b4_array != 65535) & (m_array == 1)
            formula = (b4_array - b2_array)/(b4_array + b2_array)
            ndsii = np.where(con, formula, 65535)
            try:
                sample = np.random.choice(ndsii[ndsii != 65535],50000)
                classes = jenkspy.jenks_breaks(sample, 2)
                breaks = classes[1]
                # Conversion to sea surface temperature
                formula = 1.01342+1.04948*(b20_array/100-273.15)
                b20_array = np.where(b20_array != 65535, formula, 65535)
                # Exclusion of nodata areas
                con = (b4_array != 65535) & (b20_array != 65535)
                data = con & (m_array == 1) & (ndsii != 65535)
                # Classification outcomes
                test = (b20_array <= 1) & (b4_array >= 1700) # SST + green refl.
                case1 = test & (ndsii < breaks)
                case2 = test & (ndsii >= breaks)
                case3 = (test == False) & (ndsii < breaks)
                # Classification
                # 0 = water, 1 = ice, 255 = nodata
                con = case3 | case2
                maparray = np.where(con, 255, 0)
                con = (maparray == 0) & case1
                maparray = np.where(con, 1, maparray)
                maparray = np.where(data, maparray, 255)
                if maparray.min() < 255:
                    closeraster(maparray, mod35map, ncols, nrows, gdal.GDT_Byte, 
                                255, trans)
                    print("MOD35 map created.")
                else:
                    print("MOD35 map is all nodata.")
            except:
                print("MOD35 map is all nodata.")
                
  
#****************************STEP5: MOD35 CORRECTION**************************# 
            
def step5(y,m,d):
    infolder = maskedmaps+"/"+y+"/"+m
    maplist = glob.glob(infolder+"/A"+y+d+".????.mod35map.tif")
    
    for mod35map in maplist:
        band2 = bands+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".B2.tif"
        band4 = bands+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".B4.tif"
        band7 = bands+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".B7.tif"
        band20 = bands+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".B20.tif"
        mod35 = masks+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".mod35.tif"
        correction = infolder+"/"+mod35map[-26:-13]+".mod35map_c.tif"
        tmp = workfolder+"/temp.tif"
        tmp2 = workfolder+"/temp2.tif"
        m35_array, nrows, ncols, trans, _ = openraster(mod35map, np.uint16)
        
        if not os.path.exists(correction) and np.any(m35_array == 1):
            n=0
            # Deletion of small sea ice clusters
            wbt.clump(mod35map, tmp, diag=True, zero_back=True)
            clump, nrows, ncols, trans, _ = openraster(tmp, np.int32)
            values, counts = np.unique(clump, return_counts=True)
            val_dict = dict(zip(values, counts))
            for i in val_dict.keys():
                if i < 1:
                    val_dict[i] = 65535
                elif i == 0:
                    val_dict[i] = 0
                else: 
                    if val_dict[i] >= 100:
                        val_dict[i] = 1
                    elif val_dict[i] < 100:
                        val_dict[i] = 0             
            # Map reclassification
            for i in range(nrows):
                for e in range(ncols):
                    clump[i][e] = val_dict[clump[i][e]]
            con = (clump == 0) | (clump == 65535)
            reclass = np.where(con, 0, 1)
            closeraster(reclass, tmp, ncols, nrows, gdal.GDT_UInt16, 65535, trans)
            # Calculation of distance buffer
            wbt.buffer_raster(tmp, tmp2, 35000, gridcells=False)
            b2_array, _, _, _, _ = openraster(band2, np.single)
            b4_array, _, _, _, _ = openraster(band4, np.single)
            b7_array, _, _, _, _ = openraster(band7, np.single)
            b20_array, _, _, _, _ = openraster(band20, np.single)
            tmp_array, _, _, _, _ = openraster(tmp, np.single)
            tmp2_array, _, _, _, _ = openraster(tmp2, np.single)
            m_array, _, _, _, _ = openraster(mod35, np.uint16)
            _ = None
            # Conversion of B20 to sea surface temperature
            formula = 1.01342+1.04948*(b20_array/100-273.15)
            b20_array = np.where(b20_array != 65535, formula, 65535)
            # Exclusion of nodata areas
            con = (b4_array != 65535) & (b7_array != 65535) & (b20_array != 65535)
            data = con & (m_array < 255)
            # Calculation of NDSII-2 and Jenks breaks
            formula = (b4_array - b2_array)/(b4_array + b2_array)
            ndsii = np.where(con, formula, 65535)
            # Restriction of ndsii sample to obtain Jenks breaks
            try:
                con1 = con & (b7_array <= 350) & (tmp2_array == 1)
                ndsii_s = np.where(con1, ndsii, 65535)
                sample = np.random.choice(ndsii_s[ndsii_s != 65535],50000)
                classes = jenkspy.jenks_breaks(sample, 2)
                breaks = classes[1]
            except:
                try:
                    con1 = con & (tmp2_array == 1)
                    ndsii_s = np.where(con1, ndsii, 65535)
                    sample = np.random.choice(ndsii_s[ndsii_s != 65535],50000)
                    classes = jenkspy.jenks_breaks(sample, 2)
                    breaks = classes[1]
                except:
                    # Rasters with no pixels available for correction
                    closeraster(m35_array, correction, ncols, nrows, 
                                gdal.GDT_Byte, 255, trans)
                    os.remove(tmp)
                    os.remove(tmp2)
                    print("MOD35 artefacts corrected.")
                    n=1    
            # Classification outcomes
            if n == 0:
                test = (b7_array <= 350) & (b4_array >= 1700) & (b20_array <= 1)
                test = test & (tmp2_array == 1) & (tmp_array == 0)
                test = test & (ndsii < breaks) & (m35_array != 0)
                # Classification
                # 0 = water, 1 = ice, 255 = nodata
                maparray = np.where(test, 1, 255)
                con = (maparray == 255) & (tmp_array == 1)
                maparray = np.where(con, 1, maparray)
                con = (maparray == 255) & (m35_array == 0)
                maparray = np.where(con, 0, maparray)
                maparray = np.where(data, maparray, 255)
                closeraster(maparray, correction, ncols, nrows, gdal.GDT_Byte,
                            255, trans)
                print("MOD35 artefacts corrected.")
                os.remove(tmp)
                os.remove(tmp2)
                
        elif not os.path.exists(correction) and not np.any(m35_array == 1):
            closeraster(m35_array, correction, ncols, nrows, gdal.GDT_Byte, 255, 
                        trans)
            print("MOD35 artefacts corrected.")
            
    
    
#***************************STEP6: VIS CLASSIFICATION*************************# 
            
def step6(y,m,d):
    infolder = maskedmaps+"/"+y+"/"+m
    maplist = glob.glob(infolder+"/A"+y+d+".????.mod35map.tif")
    
    for mod35map in maplist:
        band2 = bands+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".B2.tif"
        band4 = bands+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".B4.tif"
        vismask = masks+"/"+y+"/"+m+"/"+mod35map[-26:-13]+".vis.tif"
        vismap = infolder+"/"+mod35map[-26:-13]+".vismap.tif"
        
        if not os.path.exists(vismap):
            b2_array, nrows, ncols, trans, _ = openraster(band2, np.single)
            b4_array, _, _, _, _ = openraster(band4, np.single)
            v_array, _, _, _, _ = openraster(vismask, np.uint16)
            _ = None
            # Calculation of NDSII-2 and Jenks breaks
            con = (b2_array != 65535) & (b4_array != 65535) & (v_array == 1)
            formula = (b4_array - b2_array)/(b4_array + b2_array)
            ndsii = np.where(con, formula, 65535)
            try:
                sample = np.random.choice(ndsii[ndsii != 65535],50000)
                classes = jenkspy.jenks_breaks(sample, 2)
                breaks = classes[1]
                # Exclusion of nodata areas
                data = (b4_array != 65535) & (ndsii != 65535) & (v_array == 1)
                # Classification outcomes
                case1 = (b4_array >= 1700) & (ndsii < breaks)
                case2 = (b4_array >= 1700) & (ndsii >= breaks)
                case3 = ((b4_array >= 1700) == False) & (ndsii < breaks)
                # Classification
                # 0 = water, 1 = ice, 255 = nodata
                con = case1 | case2
                maparray = np.where(con, 1, 0)
                con1 = (maparray == 0) & case3
                maparray = np.where(con1, 255, maparray)
                maparray = np.where(data, maparray, 255)
                if maparray.min() < 255:
                    closeraster(maparray, vismap, ncols, nrows, gdal.GDT_Byte, 
                                255, trans)
                    print("VIS map created.")
                else:
                    print("VIS map is all nodata.")
            except:
                print("VIS map is all nodata.")
            
            
#*****************************STEP7: MAP COMPOSITING**************************#           
            
def step7(y,m,d):
    infolder = maskedmaps+"/"+y+"/"+m
    modmaplist = glob.glob(infolder+"/A"+y+d+".????.mod35map_c.tif")
    
    for mod35map in modmaplist:
        mask = masks+"/"+y+"/"+m+"/"+mod35map[-28:-15]+".mod35.tif"
        vismap = infolder+"/"+mod35map[-28:-15]+".vismap.tif"
        composite = composites+"/"+y+"/"+m+"/"+mod35map[-28:-15]+".composite.tif"

        if not os.path.exists(composite) and os.path.exists(vismap):
            m_array, nrows, ncols, trans, _ = openraster(mod35map, np.uint16)
            v_array, _, _, _, _ = openraster(vismap, np.uint16)
            m35_array, _, _, _, _ = openraster(mask, np.uint16)
            _ = None
            # Classification
            # 0 = water, 1 = ice, 254 = nodata, 255 = land
            case1 = (m_array == 1) & (v_array == 1)
            test1 = (m_array == 1) & (v_array == 0)
            test2 = (m_array == 0) & (v_array == 0)
            test3 = (m_array == 255) & (v_array == 0)
            case2 = test1 | test2 | test3
            comp_array = np.where(case1, 1, 255)
            con = (comp_array == 255) & case2
            comp_array = np.where(con, 0, comp_array)
            con = (comp_array == 255) & (m35_array < 255)
            comp_array = np.where(con, 254, comp_array)
            closeraster(comp_array, composite, ncols, nrows, gdal.GDT_Byte, 255, 
                        trans)
            print("Composite map created.")


#******************************STEP8: DAILY MAP*******************************# 
 
def step8(y,m,d):
    infolder = composites+"/"+y+"/"+m
    maplist = glob.glob(infolder+"/A"+y+d+".????.composite.tif")
    dailymap = composites+"/"+y+"/"+m+"/A"+y+d+".daily.tif"
    
    # Set the same spatial extent to all composite maps
    for i in range(len(maplist)):
        tmp = workfolder+"/tmp_"+str(i)+".tif"
        clip(maplist[i], tmp, studyarea, 255, gdal.GDT_Byte)
    maplist = glob.glob(workfolder+"/tmp_*")
    
    # Daily map creation
    if len(maplist) > 1:
        m_array, nrows, ncols, trans, _ = openraster(maplist[0], np.uint16)
        for e in maplist[1:]:
            array, _, _, _, _ = openraster(e, np.uint16)
            _ = None
            case1 = ((m_array >= 254) | (m_array == 0)) & (array == 0)
            case2 = (m_array >= 254) & (array == 254)
            case3 = ((m_array == 0) | (m_array == 1)) & (array == 1)
            case4 = (m_array >= 254) & (array == 1)
            m_array = np.where(case1, 0, m_array)
            m_array = np.where(case2, 254, m_array)
            m_array = np.where(case3, m_array+1, m_array)
            m_array = np.where(case4, 1, m_array)
        closeraster(m_array, dailymap, ncols, nrows, gdal.GDT_Byte, 255, trans)
    elif len(maplist) == 1:
        m_array, nrows, ncols, trans, _ = openraster(maplist[0], np.uint16)
        _ = None
        closeraster(m_array, dailymap, ncols, nrows, gdal.GDT_Byte, 255, trans)
        
    for e in maplist:
        os.remove(e)
    print("Daily map created.")
    
    
#*****************************STEP9: MONTHLY MAP******************************#

def step9(y,m):
    infolder = composites+"/"+y+"/"+m
    maplist = glob.glob(infolder+"/A"+y+"*daily.tif")
    likelihood = composites+"/"+y+"/"+m+"/A"+y+".likelihood_"+m+".tif"
    monthlymap = composites+"/"+y+"/"+m+"/A"+y+".month_"+m+".tif"
    tmp = workfolder+"/tmp.tif"
    
    # Creation of sea ice presence likelihood map
    if len(maplist) > 1 and not os.path.exists(likelihood):
        m_array, nrows, ncols, trans, _ = openraster(maplist[0], np.uint16)
        m_array = np.where(m_array == 254, 65534, m_array)
        m_array = np.where(m_array == 255, 65535, m_array)
        for e in maplist[1:]:
            array, _, _, _, _ = openraster(e, np.uint16)
            _ = None
            array = np.where(array == 254, 65534, array)
            array = np.where(array == 255, 65535, array)
            case1 = ((m_array >= 65534) | (m_array == 0)) & (array == 0)
            case2 = (m_array >= 65534) & (array == 65534)
            con1 = (m_array >= 0) & (m_array < 65534)
            con2 = (array >= 0) & (array < 65534)
            case3 =  con1 & con2
            case4 = (m_array >= 65534) & (array == 1)
            m_array = np.where(case1, 0, m_array)
            m_array = np.where(case2, 65534, m_array)
            m_array = np.where(case3, m_array+array, m_array)
            m_array = np.where(case4, 1, m_array)
        maximum = m_array[m_array < 65534].max()
        perc = np.round((m_array/maximum)*100, decimals=0)
        con = m_array < 65534
        m_array = np.where(con, perc, m_array)
        closeraster(m_array, likelihood, ncols, nrows, gdal.GDT_UInt16, 65535, 
                    trans)

    # Creation of monthly extent map
    if not os.path.exists(monthlymap) and len(maplist) > 1:
        m_array, nrows, ncols, trans, _ = openraster(likelihood, np.uint16)
        _ = None
        # Percent observations threshold
        con = (m_array < 10) & (m_array != 0)
        m_array = np.where(con, 65534, m_array)
        con = (m_array < 65534) & (m_array != 0)
        m_array = np.where(con, 1, m_array)
        # Reclassification
        # 0 = Nodata, 1 = Water, 2 = ice
        m_array = np.where(m_array == 1, 2, m_array)
        m_array = np.where(m_array == 0, 1, m_array)
        m_array = np.where(m_array == 65534, 0, m_array)
        closeraster(m_array, tmp, ncols, nrows, gdal.GDT_UInt16, 65535, trans)
        # Extent generation
        wbt.euclidean_allocation(tmp, monthlymap)
        os.remove(tmp)
        m_array, nrows, _, _, _ = openraster(monthlymap, np.uint16)
        _ = None
        m_array = np.where(m_array == 65535, 255, m_array)
        closeraster(m_array, monthlymap, ncols, nrows, gdal.GDT_Byte, 255, trans)
        print("Monthly map created.")

    

#***************************SCRIPT INITIALIZATION*****************************#

def run():
    directories(year_list, month_list)
    for d in days_range:
        d = str("{:03d}".format(d))
        dailymap = composites+"/"+y+"/"+m+"/A"+y+d+".daily.tif"
        if not os.path.exists(dailymap):
            step1(y,m,d)
            step2(y,m,d)
            step3(y,m,d)
            step4(y,m,d)
            step5(y,m,d)
            step6(y,m,d)
            step7(y,m,d)
            step8(y,m,d)
    monthlymap = composites+"/"+y+"/"+m+"/A"+y+".month_"+m+".tif"
    if not os.path.exists(monthlymap):
        step9(y,m)

month_days = {1:range(1,32),2:range(32,61),3:range(60,92),4:range(91,122),
              5:range(121,153),6:range(152,183),7:range(182,214),
              8:range(213,245),9:range(244,275),10:range(274,306),
              11:range(305,336),12:range(335,367)}

for y in year_list:
    for month in month_list:
        days_range = month_days[month]
        m = str("{:02d}".format(month))
        y = str(y)
        run()
        
