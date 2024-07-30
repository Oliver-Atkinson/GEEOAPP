#!/usr/bin/env python
"""Loads an image for a given region and time, runs analyses as requested and saves images to Google Drive"""

##############################################################################################################################
# Expansions:
#   Change detection - likely involved
#   Maybe find a way to choose the best (i.e. most recently launched/best resolution/most data) satellite with available data 
#   More satellites
#   Dicts of analyses types, e.g. fire detection, vegetation health, water stuff, etc.?
#   Maybe build up some kind of masking library? Something like a default snow mask that will always apply some threshold on NDSI?
##############################################################################################################################


#################################

#Introductory prints
from colorama import Style, Fore, Back
print(f"""
+------------------------------------------------+
|                                                |
|{Fore.CYAN}                  .:ttttt:.                     {Style.RESET_ALL}|
|{Fore.CYAN}                 :88S%t;t%88S:                  {Style.RESET_ALL}|
|{Fore.CYAN}               %@%.         :S8:                {Style.RESET_ALL}|
|{Fore.CYAN}             ;8;               t8.              {Style.RESET_ALL}|
|{Fore.CYAN}            :8                  .8.             {Style.RESET_ALL}|
|{Fore.CYAN}            8:                   %8             {Style.RESET_ALL}|
|{Fore.CYAN}           %S               {Fore.GREEN}..{Fore.CYAN}    @t            {Style.RESET_ALL}| 
|{Fore.CYAN}           X%               {Fore.GREEN}8X{Fore.CYAN}    %X            {Style.RESET_ALL}|
|{Fore.CYAN}           Xt              {Fore.GREEN}8tX{Fore.CYAN}    S8;           {Style.RESET_ALL}| 
|{Fore.CYAN}           t8             {Fore.GREEN}%t SX{Fore.CYAN}   8t            {Style.RESET_ALL}|
|{Fore.CYAN}            8:           {Fore.GREEN}.X. .%@{Fore.CYAN}  @             {Style.RESET_ALL}|
|{Fore.CYAN}            tS         {Fore.GREEN}.S@:   :%:{Fore.CYAN}':             {Style.RESET_ALL}|
|{Fore.CYAN}              @%       {Fore.GREEN}tX{Fore.CYAN}      ;{Fore.GREEN}t8{Fore.CYAN}              {Style.RESET_ALL}|     
|{Fore.CYAN}               S8;    {Fore.GREEN}S8{Fore.CYAN}     :%t.{Fore.GREEN}:S{Fore.CYAN}             {Style.RESET_ALL}|      
|{Fore.CYAN}                ;%XS%{Fore.GREEN}8@t{Fore.CYAN}.:::;:.   {Fore.GREEN}88:           {Style.RESET_ALL}|    
|{Fore.GREEN}                    @S;.      .;8SS8X;          {Style.RESET_ALL}|
|{Fore.GREEN}                   S8..     X8S8X:              {Style.RESET_ALL}|
|{Fore.GREEN}                  %@@8    888.::                {Style.RESET_ALL}|
|{Fore.GREEN}                :@8S8   8:;;.                   {Style.RESET_ALL}|
|{Fore.GREEN}               .:@X88Stt:;                      {Style.RESET_ALL}|
|{Fore.GREEN}              .;S88Xt;.                         {Style.RESET_ALL}|
|{Fore.GREEN}             .;88Xt;                            {Style.RESET_ALL}|
|                                                |
+------------------------------------------------+""")

print("""\n------------------------------------------------------------------------------\n
\033[1m   Welcome to the Google Earth Engine Omanos Analytics Processing Package!\033[0m
\n------------------------------------------------------------------------------\n""")

#################################

#Setting up for GEE
print("Initialising Google Earth Engine...")
import ee
ee.Initialize()
print('Earth engine initialised.')

#################################

#Other imports
import sys
from functools import partial
from IndexFunctions import *
from QuantFunctions import *
from VisFunctions import *
from CollectionLoading import *
from RunCardReader import *

logger = logging.getLogger(__name__)

print('All packages loaded.')


#################################

# Handy global variables, mainly libraries of valid processes

#Lists of indices where multiple are calculated at once
#Have to expand these every time a new analysis is added - a little clunky
NDIs_lib = ['NDVI','NDWI','NDMI','NDSI','NDRE','GNDVI']#,'NBR' - no NBR vis yet, but it can be calculated 
Ratm1s_lib = ['RECI','GCI']
SAVIs_lib = ['SAVI','OSAVI','MSAVI']
Misc_lib = ['SIPI','EVI','VARI','ARVI', 'ARI', 'MARI', 'PSRI', 'OSI', 'BSI']
Index_lib =  NDIs_lib+Ratm1s_lib+SAVIs_lib+Misc_lib+['SAR_trial']
SAR_RGB_lib = ['SAR_RGB'] 
RGB_lib = ['RGB', 'RedVegetation1','RedVegetation2', 'SWIRComposite1', 'Agriculture', 'SWIRComposite2','Bathymetric1','Bathymetric2','Bathymetric3','Bathymetric4','Bathymetric5','Bathymetric6','Urban']
RGBCalc_lib = ['Neon']
Analysis_lib = Index_lib+RGB_lib+RGBCalc_lib+SAR_RGB_lib

#List of quantative analyses
Quant_lib = ['AreaMean', 'AreaStdDev', 'AreaThreshold', 'PointValue']

#List of available operations for reducing image collections to single images
Reducers = ['median', 'mosaic']

#Allowed time units in Earth Engine (ignoring hour, minute and second)
time_units = ['year','Year', 'month','Month', 'week','Week', 'day','Day']
plural_time = [unit+'s' for unit in time_units]
time_units += plural_time

############################################

#Annoyingly, taking a mosaic or a median of a collection loses useful information that we would like to keep
#So, defining a couple of quick functions to rectify this

def mosaicer(col: ee.ImageCollection):
    """Takes a mosaic over a collection, but retains important properties, such as satellite and image dates,
    from the original collection."""
    #Get the original properties
    props = col.getInfo()['properties']
    #Take a mosaic using the in-built Earth Engine method, and imprint the properties on it
    col_mos = col.mosaic().set({'Satellite': props['Satellite'],
                'StartDate': props['StartDate'],
                'EndDate': props['EndDate'],
                'Resolution': props['Resolution']})
    return col_mos

def medianer(col: ee.ImageCollection):
    """Takes a median over a collection, but retains important properties, such as satellite and image dates,
    from the original collection."""
    #Get the original properties
    props = col.getInfo()['properties']
    #Take a median using the in-built Earth Engine method, and imprint the properties on it
    col_mos = col.median().set({'Satellite': props['Satellite'],
                'StartDate': props['StartDate'],
                'EndDate': props['EndDate'],
                'Resolution': props['Resolution']})
    return col_mos

def im_saver(im: ee.Image, analysis: str, aoi: ee.Geometry, vis: dict):
    """Saves an image to a Google Drive file.
    
    First, the filename is set, using the properties of the image to be saved and the bands chosen, with the format of
    Satellite_Bands_StartDate_EndDate. After exception handling to check that the inputs are valid, a visualisation of the
    chosen bands is made, and passed to the Earth Engine batch export method to save to Drive.
    
    Parameters:
    im: ee.Image
        The image from which to extract bands to create a visualisation, which is then saved to file
    analysis: str
        The name of the analysis which will be used to create an image - can either be a single band or three
        bands, in which case they will be saved as RGB
    aoi: ee.Geometry
        The physical area, covered by the image, that should be saved to file
    vis: dict
        Dictionary containg the minimum and maximum values (floats or ints) to be used in creating the visualisation, as well
        as the colour palette to be used (list of hex strings) for single band images
    """

    #Retrieiving the relevant image information
    im_props = im.getInfo()['properties']
    #Analysis type checking
    if type(analysis) is not str and type(analysis) is not list: raise TypeError(f'{Fore.RED}\033[1mThe analysis to be saved must an analysis name (string) or list of bands.')

    #Create a bespoke filename based on the satellite, analysis and date range
    if type(analysis) == list:
        filename = ''
        for band in analysis: 
            filename += band
        filename += "_"+im_props['Satellite'][0]+im_props['Satellite'][-1]+"_"+im_props['StartDate']+"_"+im_props['EndDate']
    else: filename = analysis+"_"+im_props['Satellite'][0]+im_props['Satellite'][-1]+"_"+im_props['StartDate']+"_"+im_props['EndDate']
    filename = aoi_name+"_"+filename 
    #Could check if the file already exists?

    #Translating the analysis to a list of strings corresponding to the bands to save    
    if analysis in RGB_lib or analysis in RGBCalc_lib or analysis in SAR_RGB_lib: bands = [analysis+'_R', analysis+'_G', analysis+'_B']
    elif type(analysis) == list: bands = analysis
    else: bands = [analysis]

    #Exception handling, otherwise code may run fine but not save the image - any other ways the image can fail?
    #Too many or too few bands chosen
    if len(bands) < 1 or len(bands) > 3: raise ValueError(f'{Fore.RED}\033[1mInvalid number of bands ("+str(len(bands))+"), must chose 1-3 bands.')
    #Ensure that the image is not empty and all chosen bands are present in the image
    im_bands = im.bandNames().getInfo()
    if len(im_bands) == 0: raise ValueError(f'{Fore.RED}\033[1mAttempted to save an empty image: '+filename)
    for band in bands: 
        if band not in im_bands: raise ValueError(f'{Fore.RED}\033[1mBand '+band+" not in image. Valid bands are: "+str(im.bandNames().getInfo()))
    
    #Values from the visualisation dictionary
    vis_min, vis_max, palette = vis['vis_min'], vis['vis_max'], vis['palette']

    #Min and max for the visualisation are invalid
    if vis_min >= vis_max: raise ValueError(f'{Fore.RED}\033[1mMinimum value greater than the maximum for the visualisation (vis_min >= vis_max)')
    #Palette is uneccessarily supplied for an 2 or 3 band image, left empty or is only a single colour for a single band image
    if len(bands) != 1 and len(palette) != 0: raise ValueError(f"{Fore.RED}\033[1mPalettes can only be specified for single band images.")
    if len(bands) == 1 and len(palette) == 1: raise ValueError(f"{Fore.RED}\033[1mPalettes must have 2 or more colours.")
    if len(bands) == 1 and len(palette) == 0: logger.warning(f"{Fore.YELLOW}Warning: No palette specified for a single band ("+bands[0]+f") image. Greyscale will be used.{Style.RESET_ALL}")

    #Creating a visualisation to save to file
    #For single band images use a given palette, if not default to greyscale (palette arg can only be given for single band images)
    if len(bands) == 1: im_save = im.visualize(**{'bands': bands, 'min': vis_min, 'max': vis_max, 'palette': palette}) 
    else: im_save = im.visualize(**{'bands': bands, 'min': vis_min, 'max': vis_max}) 
    
    #Extract the resolution to be used - may be a source of unexpected issues in certain cases (bands of different res for the same sat)
    res = im_props['Resolution']

    #Annoyingly, this method does not allow for subfolders to be used
    #The task which saves the image to Drive
    task = ee.batch.Export.image.toDrive(**{
    'image': im_save,                        #Image to be exported
    'description': filename,                 #File name
    'folder': ImSaveDir,                     #Export folder
    'scale': res,                            #Spatial resolution 
    'region': aoi                            #Region to cover
    })

    #Running the task
    task.start()

    return

def AnalysisRunner(im: ee.Image, band_dict: dict, analyses: list):
    """Takes an image and runs the requested analyses on it
    
    By making use of the functions collected in IndexFunctions, calculates a given list of indices and adds each of these
    as a band to the original image.
    
    Parameters:
    im: ee.Image
        The image for which the indices are to be calculated
    band_info: dict
        The dictionary of satellite bands that translates names to band IDs (e.g. R -> B4), as given in Satellites.yml
    analyses: list
        The list of analyses to be performed, as strings

    Returns:
    ee.Image
        The original image, with an additional band for each index calculated
    """

    #Assumes checking on analysis list has been done - see the read in function for how

    #Build lists of indices to pass to functions that calculate multiple indices
    NDIs = [ind for ind in analyses if ind in NDIs_lib]
    Ratm1s = [ind for ind in analyses if ind in Ratm1s_lib]
    SAVIs = [ind for ind in analyses if ind in SAVIs_lib]
    Composites = [comp for comp in analyses if comp in RGB_lib]

    #Creating an image to add the indices bands to
    im_plus_ind = im

    #Do the requested analyses - must be added to everytime a new analysis is added - pretty clunky but works well
    if len(NDIs) != 0: im_plus_ind = NDI_calc(im_plus_ind, band_dict, NDIs)
    if len(Ratm1s) != 0: im_plus_ind = ratm1_calc(im_plus_ind, band_dict, Ratm1s)
    if len(SAVIs) != 0: im_plus_ind = SAVIs_calc(im_plus_ind, band_dict, SAVIs)
    if len(Composites) != 0: im_plus_ind = RGBComposite(im_plus_ind, band_dict, Composites)

    if 'SIPI' in analyses: im_plus_ind = SIPI_calc(im_plus_ind, band_dict)
    if 'EVI' in analyses: im_plus_ind = EVI_calc(im_plus_ind, band_dict)
    if 'ARVI' in analyses: im_plus_ind = ARVI_calc(im_plus_ind, band_dict)
    if 'VARI' in analyses: im_plus_ind = VARI_calc(im_plus_ind, band_dict)
    if 'ARI' in analyses: im_plus_ind = ARI_calc(im_plus_ind, band_dict)
    if 'MARI' in analyses: im_plus_ind = MARI_calc(im_plus_ind, band_dict)
    if 'BSI' in analyses: im_plus_ind = BSI_calc(im_plus_ind, band_dict)
    if 'PSRI' in analyses: im_plus_ind = PSRI_calc(im_plus_ind, band_dict)
    if 'OSI' in analyses: im_plus_ind = OSI_calc(im_plus_ind, band_dict)
    if 'Neon' in analyses: im_plus_ind = Neon(im_plus_ind, band_dict)
    if 'SAR_RGB' in analyses: im_plus_ind = SAR_RGB(im_plus_ind, band_dict)
    if 'SAR_trial' in analyses: im_plus_ind = SAR_trial(im_plus_ind, band_dict)
    
    return im_plus_ind

# Maybe just return the mask so that it can then be applied to multiple images?
def ImageMasker(aoi: ee.Geometry, im: ee.Image, mask_bandName: str, upper_bound: float = 1e5, lower_bound: float = -1e5, save_raster=0, save_vector=0):
    """Applies a mask to an image from a upper and/or lower bounds on values in a given band
    
    By selecting the given band of an image, ideally an index, and setting upper and lower bounds on this index's value
    across the image, a mask is generated and applied, such that the output image is only shown for pixels that satisfy the
    bounding criteria. The mask itself can also optionally be saved as a seperate image.
    
    Parameters:
    aoi: ee.Geometry
        The physical area, covered by the image
    im: ee.Image
        The image from which to calculate, and then apply, a mask
    mask_bandName: str
        The name of the band which is used to calculate the mask
    upper_bound: float (Optional)
        The upper bound for the values in mask_bandName; values greater than this will be masked. N.B. An upper or lower bound must be passed
    lower_bound: float (Optional)
        The lower bound for the values in mask_bandName; values less than this will be masked. N.B. An upper or lower bound must be passed
    save_raster: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on saves the mask itself as a raster image (.tif)
    save_vector: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on saves the mask itself as a vector image (.shp)

    Returns:
    ee.Image
        The original image, with the mask applied to all bands."""
    
    #Exception handling takes place in the read in function

    #Getting some of the image properties 
    im_bands = im.bandNames().getInfo()
    im_props = im.getInfo()['properties']
    res = im_props['Resolution']

    #If the requested mask band is not in the image, check to see if the analysis can be calculated, if not raise an error
    if mask_bandName not in im_bands:
        if mask_bandName in Index_lib:
            im = AnalysisRunner(im, Satellite_Information[im_props['Satellite']]['Bands'], [mask_bandName])
        else: raise ValueError(f'{Fore.RED}\033[1mRequested band for mask generation, '+mask_bandName+' is neither in the image ('+str(im_bands)+') nor a valid index ('+str(Index_lib)+')')

    mask_band = im.select([mask_bandName])

    #Creating the mask
    mask_upper = mask_band.lte(upper_bound)
    mask_lower = mask_band.gte(lower_bound)
    comb_mask = mask_lower.add(mask_upper).subtract(1) #The combined mask will have values of 1 or 2; subtract 1 to shift to 0 or 1 as befits a mask

    #Check that the mask doesn't kill too much of the image
    N_pix = ee.Number(im.select(mask_bandName).reduceRegion(ee.Reducer.count(),aoi,res).get(mask_bandName))
    N_gmask = ee.Number(comb_mask.select(mask_bandName).reduceRegion(ee.Reducer.sum().unweighted(),aoi,res).get(mask_bandName))
    g_perc = N_gmask.divide(N_pix).multiply(100).getInfo()
    #Maybe include a print statement here as well, just to say what proportion of the image passes masking
    #Value of the cut offs are somewhat arbitary, and may not always be an issue, but worth having it in as a warning
    if g_perc <= 5: logger.warning(f'{Fore.YELLOW}Warning: A low proportion ('+"{0:.3g}".format(g_perc)+'%) of the image from '+im_props['StartDate']+" to "+im_props['EndDate']+f' passes the masking criteria. Images may not save correctly.{Style.RESET_ALL}')
    elif g_perc >= 95: logger.warning(f'{Fore.YELLOW}Warning: A high proportion ('+"{0:.3g}".format(g_perc)+'%) of the image from '+im_props['StartDate']+" to "+im_props['EndDate']+f' passes the masking criteria; the mask may be ineffectual.{Style.RESET_ALL}')

    #Actually apply the mask
    masked_im = im.updateMask(comb_mask)

    #If the relevant flag is switched on, save the mask, with black for masked pixels and white unmasked
    if save_raster == 1: 
        mask_save = comb_mask.visualize(**{'bands': mask_bandName, 'min': 0, 'max': 1}) 

        #Essentially just ripping out the image saving part of im_saver, with the mask as the image to visualise
        #Including the mask bounds in the filename as well, by checking if they match the defaults
        if lower_bound != -1e5 and upper_bound == 1e5: bound_str = 'gte_'+str(lower_bound)
        if upper_bound != 1e5 and lower_bound == -1e5: bound_str = 'lte_'+str(upper_bound)
        if lower_bound != -1e5 and upper_bound != 1e5: bound_str = 'gte_'+str(lower_bound)+'_lte_'+str(upper_bound)
        filename = aoi_name+'_'+mask_bandName+'_'+bound_str+'_Mask'+'_'+im_props['Satellite'][0]+im_props['Satellite'][-1]+"_"+im_props['StartDate']+"_"+im_props['EndDate']

        #The task which saves the image to Drive
        task = ee.batch.Export.image.toDrive(**{
        'image': mask_save,                       #Image to be exported
        'description': filename,                  #File name 
        'folder': MaskSaveDir,                    #Export folder 
        'scale': res,                             #Spatial resolution 
        'region': aoi                             #Region to cover
        })
    
        #Running the task
        task.start()

    #Saving as a shape file (vector) if requested
    if save_vector == 1:
        #Including the mask bounds in the filename as well, by checking if they match the defaults
        if lower_bound != -1e5 and upper_bound == 1e5: bound_str = 'gte_'+str(lower_bound)
        if upper_bound != 1e5 and lower_bound == -1e5: bound_str = 'lte_'+str(upper_bound)
        if lower_bound != -1e5 and upper_bound != 1e5: bound_str = 'gte_'+str(lower_bound)+'_lte_'+str(upper_bound)
        filename = aoi_name+'_'+mask_bandName+'_'+bound_str+'_Mask'+'_'+im_props['Satellite'][0]+im_props['Satellite'][-1]+"_"+im_props['StartDate']+"_"+im_props['EndDate']

        #Pretty much just taken from Clare's code, this exports the mask as a vector image
        #It works - not sure why things have to be the way they are, but it works
        vector = comb_mask.addBands(im.select([mask_bandName])).reduceToVectors(**{
            'geometry': ee.FeatureCollection(aoi), #Geometry just has to be a feature collection
            'crs': im.projection(),
            'scale': res,
            'geometryType': 'polygon',
            'eightConnected': False,
            'labelProperty': 'zone',
            'reducer': ee.Reducer.mean(),
        })

        task = ee.batch.Export.table.toDrive(
            collection = vector,
            description = filename, 
            folder = MaskSaveDir,
            fileFormat = 'SHP'
            )
        task.start()

    return masked_im

#Mask loading function - seems to work fairly well, if not perfectly
def MaskLoader(aoi: ee.Geometry, mask_name: str, invert: int=0, overlap_len: int=2):
    """Loads a mask from a file uploaded to Google Earth Engine.
    
    The mechanics of this function are a little unclear and the invert and overlap_len values require tuning on a file by file 
    basis, without much apparent pattern of why they should take the required values. 
    
    Parameters:
    im: ee.Image
        The image on which to apply the mask
    aoi: ee.Geometry
       The area of interest, which should span the image
    mask_name: str
        The name of the mask asset as it appears in Google Earth Engine
    invert: int [0,1]
        Flag (0/1 for off/on) to signal if the loaded mask requires inversion
    overlap_len: int [1,2] 
        The length of overlapping geometries to be clipped from the mask, should be 1 or 2

    Returns:
    ee.Image
        The mask as an Earth Engine image
    """

    print('\nLoading a mask from the shapefile saved at '+mask_name+'...')

    #Load the mask file as a feature collection
    fc = ee.FeatureCollection(mask_name)

    #Turn the feature collection into an image
    mask = ee.Image.constant(1).clip(fc.geometry()).mask().Not()
    #Get the coordinates of the geometry of this image
    coord_list = mask.getInfo()['properties']['system:footprint']['coordinates']
    #The longest coordinate list has to be removed (presumably it borders the whole aoi)
    coord_list.remove(max(coord_list, key=len))

    #Picking out multiple coordinate lists to allow for overlapping areas
    #The length should be one (returning the original list) or two, depending on unknown factors in the aoi
    overlap_areas = [ee.Geometry.Polygon(coords) for coords in coord_list if len(coords) >= overlap_len]

    #Turning the coordinate lists into a feature collections
    feature_list = [ee.Feature(ee.Geometry.Polygon(coord_list[i]), {'system:index': str(i)}) for i in range(len(coord_list))]
    feature_col = ee.FeatureCollection(feature_list)

    overlap_list = [ee.Feature(overlap_areas[i], {'system:index': str(i)}) for i in range(len(overlap_areas))]
    overlap_col = ee.FeatureCollection(overlap_list)

    #Set up base images to use
    roi_image = ee.Image(0).clip(aoi)
    loaded_image = ee.Image(1).clipToCollection(feature_col)
    #Using the 'where' method to create a first mask, then clipping to the overlaps and repeating the 'where' method
    binary_image = roi_image.where(loaded_image, loaded_image)
    overlap_image = binary_image.clipToCollection(overlap_col)
    binary_overlap = roi_image.where(overlap_image, overlap_image)

    #If the inversion flag is required, do this
    if invert == 1: binary_overlap = binary_overlap.eq(0)

    #Return the resulting mask image
    return binary_overlap

#Similarly to the masking and analyses, create some common library, e.g. healthy veg from NDVI?
def QuantAnalysis(im_list: list, areas: list, bandName: str, analyses: list, upper_bound: float = 1e5, lower_bound: float = -1e5, save_flag=1):
    """Collects the functionality for quantitative analysis over region(s) or point(s), such as mean, with plotting and illustrative image saving.
    
    The image list is iterated over to create a single array that stores all the requested analyses as a time series, with a sub array for 
    each requested analysis, which should each correspond to a predefined function. The plots for each individual region/point are then saved,
    along with a combined plot if multiple regions/points are passed to the function. The regions/points can also be highlighted and saved onto an RGB
    image if the relevant flag is switched on.

    Parameters:
    im_list: list
        A list of images, in chronological order (oldest first), in which the areas sit
    areas: list
        A list of points or regions, as ee.Point or ee.Geometry respectively, on which to perform the analysis
    bandName: str
        The name of the band on which to perform the analysis
    analyses: list
        A list of strings of analysis names, each of which should relate to a pre-defined function 
    upper_bound: float (Optional)
        Upper bound for threshold analyses
    lower_bound: float (Optional)
        Lower bound for threshold analyses
    save_flag: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on saves the most recent image with the regions/points painted on it
    """

    #The images may not be passed to this function in the correct date order when there are multiple satellites at play, so redorder to be sure
    s_dates = [im.getInfo()['properties']['StartDate'] for im in im_list]
    sorted_d = sorted(s_dates, key=lambda x: datetime.datetime.strptime(x, '%Y-%M-%d'))
    order = [s_dates.index(s_date) for s_date in sorted_d]
    im_list = [im_list[place] for place in order]

    #Get the image properties out
    ims_props = [im.getInfo()['properties'] for im in im_list]

    #-----Possible Issue------
    # Assumes first image is representative
    #If the requested mask band is not in the image, check to see if the analysis can be calculated, if not raise an error
    if bandName not in im_list[0].bandNames().getInfo():
        if bandName in Index_lib: im_list = [AnalysisRunner(im_list[i], Satellite_Information[ims_props[i]['Satellite']]['Bands'], [bandName]) for i in range(len(im_list))]
        else: raise ValueError(f'{Fore.RED}\033[1mRequested band for quantitative analysis, '+bandName+' is neither in the image ('+str(im_list[0].bandNames().getInfo())+') nor a valid index ('+str(Index_lib)+')')

    #Analysis list not a list
    if type(analyses) is not list: raise TypeError(f"{Fore.RED}\033[1mAnalysis list not a list.")
    #No analyses requested 
    if len(analyses) == 0: raise ValueError(f"{Fore.RED}\033[1mNo analyses requested - empty list passed to QuantAnalysis.")
    
    #Bounds and flag checks
    if type(lower_bound) is not int and type(lower_bound) is not float: raise TypeError(f'{Fore.RED}\033[1mLower bound for the threshold analysis must be int or float.')
    if type(upper_bound) is not int and type(upper_bound) is not float: raise TypeError(f'{Fore.RED}\033[1mUpper bound for the threshold analysis must be int or float.')
    if lower_bound >= upper_bound: raise ValueError(f"{Fore.RED}\033[1mLower bound for the threshold analysis must be below the upper bound.")
    if save_flag != 0 and save_flag != 1: raise ValueError(f"{Fore.RED}\033[1mInvalid flag value for region saving in the quantitative analysis function - only 0 and 1 accepted.")

    #Setting a name for the regions to be used in the legend, assumes all points or all areas
    if areas[0].getInfo()['type'] == 'Polygon': name = 'Area'
    else: name = 'Point'

    #Only doing point or area analyses as appropriate
    good_ans = [an for an in analyses if name in an]

    #This gets used quite a bit
    n_areas = len(areas)

    #Initalising structure to hold the data
    data = np.zeros((len(good_ans), n_areas, len(im_list)))

    s_dates, e_dates = [prop['StartDate'] for prop in ims_props], [prop['EndDate'] for prop in ims_props]
    #For each image, get the resolution and band, then perform the requested analyses
    for i in range(len(im_list)):
        im = im_list[i]
        im_props = ims_props[i]
        res = im_props['Resolution']

        plotband = im.select([bandName])

        for j in range(n_areas):
            area = areas[j]
            for k in range(len(good_ans)):
                #Threshold analyses require upper and lower bounds as additional args
                if good_ans[k] == 'AreaThreshold': data[k,j,i] = AreaThreshold(plotband, area, res, bandName, upper_bound, lower_bound)
                else: data[k,j,i] = globals()[good_ans[k]](plotband, area, res, bandName)

    #Individual plots
    for l in range(n_areas):
        for m in range(len(good_ans)):  
            analysis = good_ans[m]

            #Maybe remove if the s.d. is ever something we'd want to see seperately
            if analysis == 'AreaStdDev': continue

            #Set a filename
            filename = aoi_name+'_'+analysis+str(l+1)+'_'+bandName+"_"+s_dates[0]+'-'+e_dates[-1]

            #For plotting the mean with s.d. as an error bar
            if analysis == 'AreaMean' and 'AreaStdDev' in analyses: TimePlotter(s_dates, data[m,l,:], bandName, PlotDir+"/"+filename, errs=data[analyses.index('AreaStdDev'),l,:])
            
            #If a threshold analysis is called the filename should additionally reflect the bounds
            elif analysis == 'AreaThreshold': 
                if lower_bound != -1e5 and upper_bound == 1e5: bound_str = 'gte_'+str(lower_bound)
                if upper_bound != 1e5 and lower_bound == -1e5: bound_str = 'lte_'+str(upper_bound)
                if lower_bound != -1e5 and upper_bound != 1e5: bound_str = 'gte_'+str(lower_bound)+'_lte_'+str(upper_bound)
                filename = aoi_name+'_'+analysis+str(l+1)+'_'+bandName+'_'+bound_str+"_"+s_dates[0]+'-'+e_dates[-1]
                TimePlotter(s_dates, data[m,l,:], 'Percentage', PlotDir+"/"+filename)
            #Normal plotting otherwise
            else: TimePlotter(s_dates, data[m,l,:], bandName, PlotDir+"/"+filename)

    #Combined plot if multiple regions requested
    if n_areas > 1:
        for n in range(len(good_ans)):
            analysis = good_ans[n]
            filename = aoi_name+'_Comb'+analysis+'_'+bandName+"_"+s_dates[0]+'-'+e_dates[-1]

            #Same caveats on the plotting as above
            if analysis == 'AreaStdDev': continue
            if analysis == 'AreaMean' and 'AreaStdDev' in analyses: TimePlotter(s_dates, data[n,:,:], bandName, PlotDir+"/"+filename, name=name, errs=data[analyses.index('AreaStdDev'),:,:])
            elif analysis == 'AreaThreshold': 
                filename = aoi_name+'_Comb'+analysis+'_'+bandName+'_'+bound_str+"_"+s_dates[0]+'-'+e_dates[-1]
                TimePlotter(s_dates, data[n,:,:], 'Percentage', PlotDir+"/"+filename, name=name)
            else: TimePlotter(s_dates, data[n,:,:], bandName, PlotDir+"/"+filename, name=name)

    # Painting the regions onto the most recent RGB image and saving this, if requested
    # Would be better if the points were clearer to see, playing with colours etc.
    if save_flag == 1:
        implusvis = im
        im_sat = im.getInfo()['properties']['Satellite']
        #Slightly different painting procedures for points and areas
        if name == 'Point':
            for i in range(n_areas):
                implusvis = implusvis.paint(ee.FeatureCollection(areas[i].buffer(20)), ee.Number(0.5), ee.Number(20)) #Play with changing the width here
        elif name == 'Area':
            for i in range(n_areas):
                implusvis = implusvis.paint(ee.FeatureCollection(areas[i]), ee.Number(1), ee.Number(10))
        implusvis = implusvis.set({'Satellite': im_sat,
                    'StartDate': s_dates[-1].format('YYYY-MM-dd'),
                    'EndDate': e_dates[-1].format('YYYY-MM-dd'),
                    'Resolution': res})
        #If RGB not in the image, change this
        if 'RGB_R' not in implusvis.bandNames().getInfo(): implusvis = AnalysisRunner(implusvis, Satellite_Information[im_sat]['Bands'], ['RGB'])
        im_saver(implusvis, 'RGB', aoi, VisLoader('RGB_default', 'RGB'))

    return

def Analyser(aoi: ee.Geometry, sat: str, process_dict: dict, dates_dict: dict, mask_dict: dict, analyses_dict: dict, vis_dict: dict = {}):
    """Performs the requested analyses; getting collections and saving images along the way.
    
    Firstly, the global area of interest and date checks are performed. Runs one of the time collection generators, and creates
    a dictionary of partial functions, one for each satellite present in the collections. Iterates over the list of collections,
    reduces each of these to a single image, then performs the requested analysis and saving.

    Parameters:
    aoi: ee.Geometry
        The area of interest 
    sat: str
        The name of the satellite to use data from
    process_dict: dict
        A dictionary with the image processing parameters
            reducer: str - The name of the operation (currently median or mosaic) used to reduce collections to a single image
            cloud_cover: float - The value by which to filter for images with less than the given cloud coverage percentage
            cloud_mask: int [0,1] - Flag for application of the cloud mask
    dates_dict: dict
        A dictionary containing all the relevant date information:
            start_date: ee.Date - Start of the overall period
            end_date: ee.Date - End of the overall period
            delta_period: float (Optional) - Number of time units to step between sub-periods, e.g. sample every delta_period days/months/years, defaults to 1  
            period: str (Optional)  - The time unit to step between sub-periods, e.g. every n period, defaults to year
            delta_sub: float - Length of sub-period time to sample, in units of subperiod, e.g. a duration of delta_sub days/months/years  
            subperiod: str - Time unit for the sub-period, e.g. delta_sub subperiods
            sub_start: ee.Date - Start of the sub-period to sample (must be within the overall start and end dates), defaults to the overall start_date
    mask_dict: dict
        A dictionary containing all the relevant masking information: 
            mask_bandName: string - The name of the band which is used to calculate the mask
            upper_bound: float - Upper bound of the mask
            lower_bound: float - Lower bound of the mask
            save_raster: int - Flag (0/1 for off/on) - if on saves the mask itself as a raster image (.tif)
            save_vector: int - Flag (0/1 for off/on) - if on saves the mask itself as a vector image (.shp)
            load_name: string - Filepath for the Earth Engine asset to load as a mask
            load_invert: int [0,1] - Flag for inversion of a loaded mask
            load_overlap: int [1,2] - Length of overlap regions to clip from a loaded mask
    analyses_dict: dict
        A dictionary of the information for quantative analyses:
            reducer: string - The name of the operation to use to reduce individual collections to single images
            analyses: list - The list of indices to be calculated, as strings - includes RGB composites
            areas: list - A list of Earth Engine geometries to analyse
            bandName: string - The name of the band to perform the analysis on
            quant_analyses: list - A list of strings of analysis names, each of which should relate to a pre-defined function
            upper_bound: float - Upper bound for threshold analyses
            lower_bound: float - Lower bound for threshold analyses
            save_flag: int [0,1] - Flag (0/1 for off/on) - if on saves the most recent image with the regions/points painted on it
    vis_dict: dict
        Dictionary of visualisation dictionaries; format of {index: {vis_min: value, vis_max: value, palette: string or list of hex strings}}
    """

    #Get a list of [Collection, Satellite] pairings, using the existence of the delta sub period as a de facto flag
    if dates_dict['delta_sub'] is not None: Col_list = SubperiodSeries(aoi, sat, process_dict, dates_dict['start_date'], dates_dict['end_date'], dates_dict['delta_sub'], dates_dict['subperiod'], dates_dict['sub_start'], dates_dict['delta_period'], dates_dict['period'])
    else: Col_list = TimeSeries(aoi, sat, process_dict, dates_dict['start_date'], dates_dict['end_date'], dates_dict['delta_period'], dates_dict['period'])
    
    #Printing the image information
    n_ims = sum([ColPair[0].size().getInfo() for ColPair in Col_list])
    if len(Col_list) == 1: print('Found a total of '+'{:1g}'.format(n_ims)+' images across a single collection.')
    else: print('Found a total of '+'{:1g}'.format(n_ims)+' images across '+'{:1g}'.format(len(Col_list))+' collections.')

    #This is a bit awkward, but not horrendous - while the analysis functions can be used to map over collections, the im_saver cannot (client side requirements)
    
    #Splitting out the collection and satellite lists and reducing collections to single images
    red_list = [globals()[process_dict['reducer']](ColPair[0]) for ColPair in Col_list] 
    sat_list = [ColPair[1] for ColPair in Col_list]
    #Getting useful variables from these
    unique_sats = list(set(sat_list))
    n_ims, n_sats = len(red_list), len(unique_sats)

    #Masking options
    mask_band, mask_upper, mask_lower, mask_raster, mask_vector = mask_dict['mask_bandName'], mask_dict['upper_bound'], mask_dict['lower_bound'], mask_dict['save_raster'], mask_dict['save_vector']
    mask_file, mask_invert, mask_overlap = mask_dict['load_name'], mask_dict['load_invert'], mask_dict['load_overlap']
    if type(mask_file) is str: mask = MaskLoader(aoi, mask_file, mask_invert, mask_overlap)

    #Creating a dictionaries of partial functions, one for each satellite
    AnalysisPartials = [partial(AnalysisRunner, band_dict=Satellite_Information[u_s]['Bands'], analyses=analyses_dict['indices']) for u_s in unique_sats]
    PartialDict = dict(zip(unique_sats, AnalysisPartials))
    #Creating a dictionary of image saving functions, one for each index
    im_savers = {}
    for analysis in analyses_dict['indices']:
        #Getting the visualisation - if none specified, look up the default settings
        if analysis not in list(vis_dict.keys()): vis = Visualisations[analysis+'_default']
        else: vis = vis_dict[analysis]
        #Add the partial im_saver to the dictionary
        im_savers.update({analysis: partial(im_saver, analysis=analysis, aoi=aoi, vis=VisLoader(vis,analysis))})

    #Masking prints
    if type(mask_band) is str:
        if type(mask_file) is str: print('Generating a mask for each image that allows only pixels with values of '+mask_band+' between '+"{0:.3g}".format(mask_lower)+' and '+"{0:.3g}".format(mask_upper)+'...')
        else: print('\nGenerating a mask for each image that allows only pixels with values of '+mask_band+' between '+"{0:.3g}".format(mask_lower)+' and '+"{0:.3g}".format(mask_upper)+'...')
    
    im_flow = [] #List to carry images through the flow
    #Do the required masking, if any
    for im in red_list:
        if type(mask_band) is str: im = ImageMasker(aoi, im, mask_band, mask_upper, mask_lower, mask_raster, mask_vector)
        if type(mask_file) is str: im = im.updateMask(mask)
        im_flow.append(im)
    if type(mask_band) is str or type(mask_file) is str: print("Masking completed.")

    #For each satellite, use the list of masked and reduced images for the analyses with the partial functions
    im_flow2 = [] #A second list to carry images through the flow
    if len(AnalysisPartials) != 0: print("\nCalculating the requested indices...")
    for i in range(n_sats):
        col_sat = unique_sats[i]
        Analysed_list = [PartialDict[col_sat](im_flow[j]) for j in range(n_ims) if sat_list[j] == col_sat]
        #Iterate over each image in the satellite list, saving each analysis requested
        for im in Analysed_list:
            for analysis in analyses_dict['indices']:
                im_savers[analysis](im)
            im_flow2.append(im) 
    if len(AnalysisPartials) != 0: print("Index calculations completed.")

    #If there are no analyses, just pass the original flow list
    if len(im_flow2) == 0: im_flow2 = im_flow

    #Do the quantative analyses on the full selection of images 
    if len(analyses_dict['quant_analyses']) != 0: 
        print('\nPerforming the requested quantitative analyses...')
        #Doing points and areas in seperate runs if both types of geometry are requested
        if len(analyses_dict['areas'][0]) > 0: QuantAnalysis(im_flow2, analyses_dict['areas'][0], analyses_dict['bandName'], analyses_dict['quant_analyses'], analyses_dict['upper_bound'], analyses_dict['lower_bound'], analyses_dict['save_flag'])
        if len(analyses_dict['areas'][1]) > 0: QuantAnalysis(im_flow2, analyses_dict['areas'][1], analyses_dict['bandName'], analyses_dict['quant_analyses'], analyses_dict['upper_bound'], analyses_dict['lower_bound'], analyses_dict['save_flag'])
        print('Quantitative analyses completed.')

    return


######################################

#Actually running everything

#################################

#Identifying the run card to use from the command line input
if len(sys.argv) == 2:
    RunCard = sys.argv[1]

    #Read all the information in 
    print('\nLoading and checking user inputs from '+RunCard+'...')
    names, aoi, sat, process_dict, dates_dict, mask_dict, analyses_dict, vis_dict = RunCardReader(RunCard)
    print(f'{Fore.GREEN}User input read in successfully from '+RunCard+f".{Style.RESET_ALL}\n")

#If no run card is given, generate one using the GUI
else:
    print('\nInitialising a GUI for user input...')
    from GUI import *
    print(f'{Fore.GREEN}GUI input validated!\n{Style.RESET_ALL}')

#Setting useful global variables
aoi_name, ImSaveDir, MaskSaveDir, PlotDir = names['Name'], names['ImSaveDir'], names['MaskSaveDir'], names['PlotDir']

#Finally, doing the analysis
Analyser(aoi, sat, process_dict, dates_dict, mask_dict, analyses_dict, vis_dict)
print(f'{Fore.GREEN}\n\033[1mAll computations completed!\033[1m\n')

######################################