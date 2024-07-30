#!/usr/bin/env python
"""Collects all the functions used for the numerous indices.
These take an image and return it with a band added for each index requested.
Where there are multiple similar indices (e.g. NDIs) these are calculated in the same function, taking a list of index names."""

#Imports
import ee
from colorama import Style, Fore, Back


###############################################################
# To do:
#   Any more indices to add
#       - https://github.com/sentinel-hub/custom-scripts/tree/master/sentinel-2/apa_script
#       - https://github.com/sentinel-hub/custom-scripts/tree/master/sentinel-2/indexdb - index index
#       - https://github.com/sentinel-hub/custom-scripts/blob/master/sentinel-2/barren_soil - also has rgb vis
#       - https://github.com/sentinel-hub/custom-scripts/tree/master/sentinel-2/oil-spill-index - also has rgb vis
#       - https://github.com/sentinel-hub/custom-scripts/blob/master/sentinel-2/urban_classified/script.js
#       - https://github.com/sentinel-hub/custom-scripts/tree/master/sentinel-2/water_surface_visualizer
#   Exception handling - any other fail cases?
#   Maybe change actual calculations to the expression method?
#   https://github.com/sentinel-hub/custom-scripts/tree/master/sentinel-2/composites - useful info on 3 band images
###############################################################


def BandChecker(sat_bands: dict, name: str, bands: list):
    """
    Checks that given bands are indeed present in the image, by iterating over
    a list of bands and checking each one is present in the satellite data, and if not, raising an error.
    
    Parameters:
    sat_bands: dict
        The dictionary of satellite bands that translates names to band IDs (e.g. R -> B4), as given in Satellites.yml
    name: str
        The name of the index being checked, only used for a more complete error message
    bands: list
        A list of strings, each corresponding to a band, the presence of which in sat_bands will be checked
    """
    
    for band in bands: 
        try: sat_bands[band]
        except KeyError: raise ValueError(f"{Fore.RED}\033[1mBand "+band+", required for "+name+", not in image.")

    return 

def NDI_calc(im: ee.Image, band_dict: dict, NDIs: list):
    """
    Calculates a list of Normalised Differences.

    Iterating over the given list of NDIs, extract the relevant bands required from the hardcoded dictionary,
    check that these are included in the satellite's bands. The NDIs are then calculated using the built-in
    Earth Engine method, and these NDI bands are added to the original image.
    
    Parameters:
    im: ee.Image
        The image for which the NDIs are to be calculated
    band_dict: dict
        The dictionary of satellite bands that translates names to band IDs (e.g. R -> B4), as given in Satellites.yml
    NDIs: list
        A list of NDIs to be calculated, as strings

    Returns:
    ee.Image
        The original image, with an additional band for each NDI calculated
    """
    
    #---------Possible issue with band names, e.g. SWIR1---------
    #First, a dictionary of NDIS by the bands they use (same ordering as the GEE function, i.e. (B1 - B2) / (B1 + B2))
    NDI_dict = {'NDVI': ['NIR', 'R'],
                'NDWI': ['G', 'NIR'],
                'NDMI': ['NIR', 'SWIR1'], 
                'NDSI': ['G', 'SWIR1'], 
                'NBR': ['R', 'SWIR2'],   
                'NDRE': ['NIR', 'RE1'], 
                'GNDVI': ['NIR', 'G']}

    #Exception handling
    #Zero length list of NDIs passed
    if len(NDIs) == 0: raise ValueError(f"{Fore.RED}\033[1mAttemping to calculate an empty list of NDIs.")
    for NDI in NDIs:
        #Check that the requested NDIs are valid
        try: NDI_dict[NDI]
        except KeyError: raise ValueError(f'{Fore.RED}\033[1m'+NDI+" not recognised, currently supported NDIS are "+str(list(NDI_dict.keys())))
        #Ensure that all chosen bands are in the image 
        BandChecker(band_dict, NDI, NDI_dict[NDI]) #Will lead to repeats when multiple NDIs are calculated

    #Iterating over the requested NDIs
    NDI_bands = [im.normalizedDifference([band_dict[NDI_dict[NDI][0]], band_dict[NDI_dict[NDI][1]]]).rename(NDI) for NDI in NDIs]

    #Adding the NDI bands to the original image
    NDI_bands_im = ee.Image(NDI_bands)
    im_plus_NDIs = im.addBands(NDI_bands_im)

    return im_plus_NDIs

def ratm1_calc(im: ee.Image, band_dict: dict, rats: list):
    """
    Calculates a list of 'ratio minus ones' ((B1/B2) - 1)).

    Iterating over the given list of ratios, extract the relevant bands required from the hardcoded dictionary,
    check that these are included in the satellite's bands. The ratio minus ones are then calculated manually
    and these bands are added to the original image.
    
    Parameters:
    im: ee.Image
        The image for which the ratio minus ones are to be calculated
    band_dict: dict
        The dictionary of satellite bands that translates names to band IDs (e.g. R -> B4), as given in Satellites.yml
    rats: list
        A list of ratio minus ones to be calculated, as strings

    Returns:
    ee.Image
        The original image, with an additional band for each ratio minus one calculated
    """
    
    #Dictionary of possible cases - probably will expand
    rat_dict = {'RECI': ['NIR', 'R'],
                'GCI': ['NIR', 'G']}

    #Exception handling
    #Zero length list of ratios passed
    if len(rats) == 0: raise ValueError(f"{Fore.RED}\033[1mAttemping to calculate an empty list of ratio minus ones.")
    for rat in rats:
        #Check that the requested ratios are valid
        try: rat_dict[rat]
        except KeyError: raise ValueError(f'{Fore.RED}\033[1m'+rat+" not recognised, currently supported ratio minus ones are "+str(list(rat_dict.keys())))
        #Ensure that all chosen bands are in the image 
        BandChecker(band_dict, rat, rat_dict[rat])

    #Doing the the ratio calculation
    #List comprehension would make this snappier but also quite a girthy line
    rat_bands = []
    for rat in rats:
        num, denom = im.select(band_dict[rat_dict[rat][0]]), im.select(band_dict[rat_dict[rat][1]])
        ratm1 = num.divide(denom).subtract(1).rename(rat)
        rat_bands.append(ratm1)

    #Making an image out of the bands and returning the original image with these bands added
    rat_im = ee.Image(rat_bands)

    return im.addBands(rat_im)
    
def SIPI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Structure Intensive Pigment Vegetation Index; (NIR - B) / (NIR - R), and adds a corresponding
    band to the image"""

    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'SIPI', ['NIR', 'R', 'B'])

    #Getting the bands
    NIR, R, B = im.select(band_dict["NIR"]), im.select(band_dict["R"]), im.select(band_dict["B"]) 

    #Doing the calculation - splitting it into parts to make it more readable
    num = NIR.subtract(B)
    denom = NIR.subtract(R)
    SIPI = num.divide(denom).rename("SIPI")

    #Adding the EVI band to the image
    return im.addBands(SIPI) 

def EVI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Enhanced Vegetation Index with a fixed set of parameters G, C1, C2, L, and adds a corresponding
    band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'EVI', ['NIR', 'R', 'B'])

    #Defining constants, using "typical" values (for now at least)
    G, C1, C2, L = 2.5, 6, 7.5, 1

    #Getting the bands
    NIR, R, B = im.select(band_dict["NIR"]), im.select(band_dict["R"]), im.select(band_dict["B"]) 

    #Doing the calculation - splitting it into parts to make it more readable
    num = NIR.subtract(R)
    denom = NIR.add(R.multiply(C1)).subtract(B.multiply(C2)).add(L)
    EVI = num.divide(denom).multiply(G).rename("EVI")

    #Adding the EVI band to the image
    return im.addBands(EVI)

def SAVIs_calc(im: ee.Image, band_dict: dict, SAVIs: list):
    """Calculates the Soil Adjusted Vegetation Indices, including the Modified and Optimised versions with fixed L = 0.5,
    and adds these as bands to the original image. The form of these indices is more involved than most, but it makes
    most sense to group all these indices together in a single function.
    
    Parameters:
    im: ee.Image
        The image for which the Soil Adjusted Vegetation Indices are to be calculated
    band_dict: dict
        The dictionary of satellite bands that translates names to band IDs (e.g. R -> B4), as given in Satellites.yml
    SAVIs: list
        A list of Soil Adjusted Vegetation Indices to be calculated, as strings

    Returns:
    ee.Image
        The original image, with an additional band for each Soil Adjusted Vegetation Index calculated
    """
    
    #Exception handling
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'SAVIs', ['NIR', 'R'])
    #Zero length list of SAVIs passed
    if len(SAVIs) == 0: raise ValueError(f"{Fore.RED}\033[1mAttemping to calculate an empty list of SAVIs.")
    for sav in SAVIs:
        #Check that the requested SAVIs are valid
        if sav not in ['SAVI', 'OSAVI', 'MSAVI']: raise ValueError(f'{Fore.RED}\033[1m'+sav+" not recognised, currently supported SAVIs are SAVI, OSAVI and MSAVI")

    #Getting the bands
    NIR, R = im.select(band_dict["NIR"]), im.select(band_dict["R"])

    #Defining L, typically 0.5, otherwise -1 <= L <= 1 
    L = 0.5

    #Creating an empty image to hold the SAVIs
    SAVI_im = ee.Image()

    #Doing the SAVI calculation, if required
    if 'SAVI' in SAVIs:
        S_num = NIR.subtract(R).multiply(1+L)
        S_denom = NIR.add(R).add(L)
        SAVI = S_num.divide(S_denom).rename("SAVI")
        SAVI_im = SAVI_im.addBands(SAVI)

    #Doing the optimised SAVI calculation, if required
    if 'OSAVI' in SAVIs:
        O_num = NIR.subtract(R)
        O_denom = NIR.add(R).add(0.16) #0.16 by definition, other values may be used
        OSAVI = O_num.divide(O_denom).rename("OSAVI")
        SAVI_im = SAVI_im.addBands(OSAVI)
    
    #Modified SAVI, going termwise: MSAVI = 0.5 * (T1 - T2), T1 = 2*NIR + 1, T2 = sqrt((2*NIR+1)^2 - 8*(NIR - R))
    if 'MSAVI' in SAVIs:
        T1 = NIR.multiply(2).add(1)
        T2_0 = T1.pow(2)
        T2_1 = NIR.subtract(R).multiply(8)
        T2 = (T2_0.subtract(T2_1)).sqrt()
        MSAVI = (T1.subtract(T2)).multiply(0.5).rename("MSAVI")
        SAVI_im = SAVI_im.addBands(MSAVI)

    #Adding the SAVI bands to the original image
    return im.addBands(SAVI_im.select(SAVIs)) #Using select to avoid the constant band the empty image constructor

def ARVI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Atmoshperically Resistant Vegetation Index (NIR - 2*R + B)/(NIR + 2*R + B), and adds a corresponding band to the image"""

    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'ARVI', ['NIR', 'R', 'B'])

    #Getting the bands
    NIR, R, B = im.select(band_dict["NIR"]), im.select(band_dict["R"]), im.select(band_dict["B"]) 

    #Doing the calculation
    num = NIR.subtract(R.multiply(2)).add(B)
    denom = NIR.add(R.multiply(2)).add(B)
    ARVI = num.divide(denom).rename("ARVI")

    #Adding the ARVI band to the image
    return im.addBands(ARVI)

def VARI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Visible Atmoshperically Resistant Vegetation Index (G - R)/(G + R - B), and adds a corresponding band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'VARI', ['R', 'G', 'B'])

    #Getting the bands
    R, G, B = im.select(band_dict["R"]), im.select(band_dict["G"]), im.select(band_dict["B"]) 

    #Doing the calculation
    num = G.subtract(R)
    denom = G.add(R).subtract(B)
    VARI = num.divide(denom).rename("VARI")

    #Adding the VARI band to the image
    return im.addBands(VARI)

def ARI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Anthocyanin Reflectance Index 1/G - 1/RE1, and adds a corresponding band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'ARI', ['G', 'RE1'])

    #Getting the bands
    G, RE1 = im.select(band_dict["G"]), im.select(band_dict["RE1"])

    #Doing the calculation
    T1 = ee.Image.constant(1).divide(G)
    T2 = ee.Image.constant(1).divide(RE1)
    ARI = T1.subtract(T2).rename("ARI")

    #Adding the ARI band to the image
    return im.addBands(ARI)

def MARI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Modified Anthocyanin Reflectance Index (1/G - 1/RE1) * NIR, and adds a corresponding band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'MARI', ['G', 'RE1', 'NIR'])

    #Getting the bands
    G, RE1, NIR = im.select(band_dict["G"]), im.select(band_dict["RE1"]), im.select(band_dict["NIR"])

    #Doing the calculation
    T1 = ee.Image.constant(1).divide(G)
    T2 = ee.Image.constant(1).divide(RE1)
    ARI = T1.subtract(T2)
    MARI = ARI.multiply(NIR).rename("MARI")

    #Adding the ARI band to the image
    return im.addBands(MARI)

def BSI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Barren Soil Index (SWIR2+R-NIR-B)/(SWIR2+R+NIR+B), and adds a corresponding band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'BSI', ['SWIR1', 'R', 'NIR', 'B'])

    #Getting the bands
    SWIR1, R, NIR, B = im.select(band_dict["SWIR1"]), im.select(band_dict["R"]), im.select(band_dict["NIR"]), im.select(band_dict["B"])

    #Doing the calculation
    num = SWIR1.add(R).subtract(NIR).subtract(B)
    denom = SWIR1.add(R).add(NIR).add(B)
    BSI = num.divide(denom).rename("BSI")

    #Adding the BSI band to the image
    return im.addBands(BSI)
    
def PSRI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Plant Senescence Reflectance Index (R-B)/RE2, and adds a corresponding band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'PSRI', ['R', 'B', 'RE2'])

    #Getting the bands
    R, B, RE2 = im.select(band_dict["R"]), im.select(band_dict["B"]), im.select(band_dict["RE2"])

    #Doing the calculation
    num = R.subtract(B)
    PSRI = num.divide(RE2).rename("PSRI")

    #Adding the BSI band to the image
    return im.addBands(PSRI)

def OSI_calc(im: ee.Image, band_dict: dict):
    """Calculates the Oil Spill Index (G+R)/B, and adds a corresponding band to the image"""
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'OSI', ['G', 'R', 'B'])

    #Getting the bands
    G, R, B = im.select(band_dict["G"]), im.select(band_dict["R"]), im.select(band_dict["B"])

    #Doing the calculation
    num = G.add(R)
    OSI = num.divide(B).rename("OSI")

    #Adding the BSI band to the image
    return im.addBands(OSI)

#---------Minor issue---------
#Doing this duplication of the bands, especially if repeated for several three band visualisations, may lead to very large/unweildy images 
def RGBComposite(im: ee.Image, band_dict: dict, comps: list):
    """"Handles the RGB composite image saving"""

    #Exception handling
    #Zero length list of NDIs passed
    if len(comps) == 0: raise ValueError(f"{Fore.RED}\033[1mAttemping to create an empty list of RGB composites.")

    #First, a dictionary of composites by the bands they use (in RGB order)
    Comp_dict = {'RGB': ['R', 'G', 'B'],                    #True colour
                'RedVegetation1': ['NIR', 'G', 'R'],        #Highlights healthy vegetation in red
                'RedVegetation2': ['NIR', 'R', 'B'],
                'SWIRComposite1': ['SWIR2', 'NIR', 'R'], 
                'Agriculture': ['SWIR1', 'NIR', 'B'], 
                'SWIRComposite2': ['SWIR2', 'SWIR1', 'B'],  
                'Bathymetric1': ['R', 'B', 'Aerosols'],
                'Bathymetric2': ['NIR', 'RE2', 'R'],
                'Bathymetric3': ['NIR', 'RE1', 'R'],
                'Bathymetric4': ['NIR', 'SWIR1', 'R'],
                'Bathymetric5': ['NIR', 'SWIR1', 'SWIR2'],
                'Bathymetric6': ['SWIR1', 'NIR', 'G'],
                'Urban': ['SWIR2', 'SWIR1', 'R']
                }

    #Holding image
    im_plusRGBs = im
    #Iterate over the requested composites and add the RGB bands to the image
    for comp in comps:
        #Check that the requested NDIs are valid
        try: Comp_dict[comp]
        except KeyError: raise ValueError(f"{Fore.RED}\033[1m"+comp+" not recognised, currently supported false colour composites are "+str(list(Comp_dict.keys())))
        #Ensure that all chosen bands are in the image 
        BandChecker(band_dict, comp, Comp_dict[comp]) #Will lead to repeats when multiple composites are made

        R_name, G_name, B_name = Comp_dict[comp][0], Comp_dict[comp][1], Comp_dict[comp][2]
        RGB = [im.select(band_dict[R_name]).rename(comp+'_R'), im.select(band_dict[G_name]).rename(comp+'_G'), im.select(band_dict[B_name]).rename(comp+'_B')]
        im_plusRGBs = im_plusRGBs.addBands(RGB)

    return im_plusRGBs

def Neon(im: ee.Image, band_dict: dict):
    
    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'Neon', ['R', 'G', 'B', 'SWIR2'])

    #Getting the bands
    R, G, B, SWIR = im.select(band_dict["R"]), im.select(band_dict["G"]), im.select(band_dict["B"]), im.select(band_dict["SWIR2"])

    # gain, gam = 2.3, -0.95 #recommended gamma: -0.55 to -0.95

    Val_R = SWIR.subtract(R).multiply(3.8).rename('Neon_R')
    Val_G = G.multiply(2.3).rename('Neon_G')  
    Val_B = B.multiply(4.2).rename('Neon_B')   

    return im.addBands([Val_R, Val_G, Val_B])

#Quick RGB composite from the S1 change detection tutorial
def SAR_RGB(im: ee.Image, band_dict: dict):

    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'SAR RGB', ['VV', 'VH'])

    R, G = im.select('VV').rename('SAR_RGB_R'), im.select('VH').rename('SAR_RGB_G')
    B = R.divide(G).rename('SAR_RGB_B')

    return im.addBands([R, G, B])

#Trying to tease out above ground biomass from SAR data
def SAR_trial(im: ee.Image, band_dict: dict):

    #Ensure that all chosen bands are in the image 
    BandChecker(band_dict, 'SAR trial', ['VV', 'VH'])

    VV, VH = im.select('VH'), im.select('VV')
    # abg = im.select('VH').divide(im.select('VV')).rename('SAR_trial')
    # abg = im.normalizedDifference(['VH', 'VV']).rename('SAR_trial')
    num, denom = VV.subtract(VH), VV.add(VH)
    abg = num.divide(denom).rename('SAR_trial')


    aoi = ee.Geometry.Rectangle([[3.3155, 36.6850], [3.6025, 36.8162]])
    min_max= abg.reduceRegion(ee.Reducer.minMax(),aoi,10).getInfo()

    return im.addBands([abg])