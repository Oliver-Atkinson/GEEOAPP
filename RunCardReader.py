"""Functionality for reading in the run vard and checking all the user inputs"""

#Inports
import ee
import datetime
import os
from VisFunctions import *

logger = logging.getLogger(__name__)


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

# Reading in the satellite information
with open("Satellites.yml", "r") as stream:
    try: Satellite_Information = yaml.safe_load(stream)
    except yaml.YAMLError as exc: print(f'{Fore.RED}\033[1m'+str(exc))

############################################

#Exception handling function for all flags - if not included, set to off
def FlagChecker(flag, name):
    if flag is None: flag = 0
    if flag != 0 and flag != 1: raise ValueError(f"{Fore.RED}\033[1mInvalid value for "+name+f" flag - only 0 and 1 accepted.{Style.RESET_ALL}")
    return flag

#Checking if numbers are indeed valid numbers (and converting to floats as required)
def NumChecker(value, name):
    if value is None: return #May be an issue if some required numbers are not included in the run card
    try: value = float(value)
    except ValueError: raise TypeError(f"{Fore.RED}\033[1m"+name+f' must be a number.{Style.RESET_ALL}')
    return value


#Checks for valid coordinates
def CoordChecker(value, name):
        NumChecker(value, name)
        if 'long' in name:
            if abs(value) > 180: raise ValueError(f'{Fore.RED}\033[1mLongitudes must be between -180 and 180 degrees; given value was '+str(value)+f'{Style.RESET_ALL}')
        if 'lat' in name:
            if abs(value) > 90: raise ValueError(f'{Fore.RED}\033[1mLatitudes must be between -90 and 90 degrees; given value was '+str(value)+f'{Style.RESET_ALL}')
        return

#Type and valid character check on the names
def NameChecker(filename, name):
    if type(filename) is not str: raise TypeError(f"{Fore.RED}\033[1m"+name+f' name must be a string.{Style.RESET_ALL}')
    safe_name = re.sub(r'[^\w_. -]', '', filename)
    if safe_name != filename: logger.warning(f'{Fore.YELLOW}Warning: Unsafe characters in '+name+' name. Renamed to '+safe_name+f'{Style.RESET_ALL}')
    return safe_name  

def DateChecker(sat: str, start_date: ee.Date, end_date: ee.Date):
    """Checks that the dates given make sense; i.e. are valid date formats, are not in the future, that the start date is
    before the end date and that the satellite from which the data is requested was operational for at least some of the
    period requested."""

    #Requested range
    req_range = ee.DateRange(start_date, end_date)
    #Checking that the dates are in fact dates
    try: start_date.getInfo()
    except ee.ee_exception.EEException: raise ValueError(f"{Fore.RED}\033[1mInvalid start date format; use 'YYYY-MM-dd'.{Style.RESET_ALL}")
    try: end_date.getInfo()
    except ee.ee_exception.EEException: raise ValueError(f"{Fore.RED}\033[1mInvalid end date format; use 'YYYY-MM-dd'.{Style.RESET_ALL}")
    #Check that the start date is before the end date
    if req_range.isEmpty().getInfo() == True: raise ValueError(f"{Fore.RED}\033[1mEmpty date range, most likely because the start date is after the end date.{Style.RESET_ALL}")
    
    #Check that the dates are not in the future (defined as from tomorrow onwards)
    today = ee.Date(str(datetime.date.today())) #Get todays date
    tomorrow = today.advance(1, "day")
    if tomorrow.difference(start_date, 'day').getInfo() < 0: raise ValueError(f"{Fore.RED}\033[1mStart date in the future.{Style.RESET_ALL}")
    if tomorrow.difference(end_date, 'day').getInfo() < 0: raise ValueError(f"{Fore.RED}\033[1mEnd date in the future.{Style.RESET_ALL}")

    #Checking that the satellite was operational during the dates requested
    sat_dates = Satellite_Information[sat]['Operative']
    sat_start, sat_end = sat_dates[0], sat_dates[1]
    if sat_dates[1] == 'Present': sat_end = today #If the satellite is still active, set the end date to today - not always accurate as data not up to the exact date
    sat_range = ee.DateRange(ee.Date(sat_start), ee.Date(sat_end))

    #Check if there's an overlap, and if both the start and end dates are included in the satellite dates
    if req_range.intersection(sat_range).isEmpty().getInfo() == True: raise ValueError(f"{Fore.RED}\033[1mNo overlap between the requested dates and "+sat+" dates of operation ("+sat_dates[0]+' - '+sat_dates[1]+f').{Style.RESET_ALL}')
    if sat_range.contains(start_date).getInfo() == False: logger.warning(f"{Fore.YELLOW}Warning: Start date is not during "+sat+" operation ("+sat_dates[0]+' - '+sat_dates[1]+f').{Style.RESET_ALL}')
    if sat_range.contains(end_date).getInfo() == False: logger.warning(f"{Fore.YELLOW}Warning: End date is not during "+sat+" operation ("+sat_dates[0]+' - '+sat_dates[1]+f').{Style.RESET_ALL}')

    #Extra issue of Landsat 4 and the data loss (see https://web.archive.org/web/20100821074100/http://landsat.gsfc.nasa.gov/about/landsat4.html)
    if sat == 'Landsat 4':
        L4_loss = ee.DateRange(ee.Date('1984-04-01'), ee.Date('1987-05-31'))
        if req_range.intersection(L4_loss).isEmpty().getInfo() == False: logger.warning(f"{Fore.YELLOW}Landsat 4 suffered data loss between April 1984 and May 1987, which overlaps with the requested dates.{Style.RESET_ALL}")

    #Issues in the Landsat 7 dataset from May 2003 onwards
    if sat == 'Landsat 7':
        L7_loss = ee.DateRange(ee.Date('2003-05-01'), ee.Date('2022-08-15'))
        if req_range.intersection(L7_loss).isEmpty().getInfo() == False: logger.warning(f"{Fore.YELLOW}Landsat 7 data from after May 2003 has issues; images may be lacking stripes of data.{Style.RESET_ALL}")

    return

def AoiChecker(sat: str, aoi: ee.Geometry):
    """Checks that the area of interest is valid, i.e. is of the right form and is not too large or too small."""

    res = Satellite_Information[sat]['Resolution']

    #Check that the aoi is formatted correctly - should be removed once user input chain is more robust/exists
    if aoi.getInfo()['type'] != 'Polygon': raise TypeError(f"{Fore.RED}\033[1mArea of interest not formatted as a valid shape.{Style.RESET_ALL}")
    #Check that the aoi is bounded
    if aoi.isUnbounded().getInfo() == True: raise ValueError(f"{Fore.RED}\033[1mArea of interest is unbounded.{Style.RESET_ALL}")

    area = aoi.area().getInfo() #Surface area in m^2
    ratio = area/(res**2) #Ratio of physical area to pixel area ~ N of pixels
    #If the area is exceedingly large, don't even try to run it - numbers are semi-arbitary for the large cases
    if ratio > 30e6: 
        logger.warning(f"{Fore.YELLOW}""Warning: The area of interest is very large compared to the satellite resolution; images will take exceedingly long times to generate, \
so this run will now be terminated. Consider refining the area of interest or choosing a satellite with lower resolution."""+f'{Style.RESET_ALL}')
        # exit()
    #Too large; files will save but will take some time, so raise a warning - somewhat arbitary threshold
    if ratio > 5e6: logger.warning(f"{Fore.YELLOW}Warning: The area of interest is large compared to the satellite resolution; images may take some time to generate.{Style.RESET_ALL}")
    #If the physical area is too small, compared to the satellite resolution, no image will be returned
    if ratio < 38e3: raise ValueError(f"{Fore.RED}\033[1mArea of interest is too small to return an image - expand or choose a higher resolution satellite.{Style.RESET_ALL}")

    return

#Open up the runcard file and extract all the relevant information, with plenty of input checking done 
def RunCardReader(Runcard: str):
    #Opening up the given run card and extracting the full data
    with open('RunCards/'+Runcard, "r") as stream:
        try: data = yaml.safe_load(stream)
        except yaml.YAMLError as exc: raise ValueError(f"{Fore.RED}\033[1m"+str(exc))
    
    #Getting the satellite out
    sat = data['Satellite']
    if type(sat) is not str: raise TypeError(f'{Fore.RED}\033[1mSatellite name must be a string.{Style.RESET_ALL}')
    if sat not in Satellite_Information.keys(): raise ValueError(f"{Fore.RED}\033[1m"+str(sat)+" not recognised, currently supported satellites are "+str(list(Satellite_Information.keys()))+f'{Style.RESET_ALL}')

    #Splitting out the dictionaries to individual ones
    names, process, dates, masking, quants, vis = data['Names'], data['Processing'], data['Dates'], data['Masking'], data['Analyses'], data['Visualisations']
    
    #Checking names
    names['Name'], names['ImSaveDir'], names['MaskSaveDir'], names['PlotDir'] = NameChecker(names['Name'], 'AOI'), NameChecker(names['ImSaveDir'], 'Image saving folder'), NameChecker(names['MaskSaveDir'], 'Mask saving folder'), NameChecker(names['PlotDir'], 'Plot saving folder')
    #Check that local directories exist, if not, make them
    if os.path.isdir(names['PlotDir']) == False: os.makedirs(names['PlotDir'])

    #Image processing 
    if process['reducer'] is None: raise ValueError(f"{Fore.RED}\033[1mRequested reducer not recognised. Valid reducers are "+str(Reducers)+f'{Style.RESET_ALL}')
    if process['reducer'].lower() not in Reducers: raise ValueError(f"{Fore.RED}\033[1mRequested reducer not recognised. Valid reducers are "+str(Reducers)+f'{Style.RESET_ALL}')
    process['reducer'] = process['reducer'].lower()+'er'
    # Cloud perc must be number between 0 and 100, if None set to 100
    if process['cloud_cover'] is None: process['cloud_cover'] = 100
    process['cloud_cover'] = NumChecker(process['cloud_cover'], 'Maximum cloud cover percentage')
    if process['cloud_cover'] > 100 or process['cloud_cover'] < 0: raise ValueError(f'{Fore.RED}\033[1mMaximum cloud cover percentage must be between 0 and 100.{Style.RESET_ALL}')
    process['cloud_mask'] = FlagChecker(process['cloud_mask'], 'Cloud mask')
    
    #Checking validity of the coordinates
    for coord in ['min_long', 'min_lat', 'max_long', 'max_lat']:
        CoordChecker(data[coord], coord)

    #Changing the aoi coordinates to an EE geometry - assume rectangular area
    aoi = ee.Geometry.Rectangle(data['min_long'], data['min_lat'], data['max_long'], data['max_lat'])
    #Global AOI checker
    AoiChecker(sat, aoi)

    #Date dictionary; converting to Earth Engine date objects and error handling
    dates['start_date'], dates['end_date'] = ee.Date(dates['start_date']), ee.Date(dates['end_date']) #Start and end date checks are handled by the DateChecker function, run in Analyser
    #Exception handling for the subperiod cases
    if dates['delta_sub'] is not None:
        if dates['sub_start'] is None: dates['sub_start'] = 0 #Set default sub start 
        if dates['subperiod'] not in time_units: raise ValueError(f'{Fore.RED}\033[1mSub-period must be a valid time unit.{Style.RESET_ALL}')
    if dates['sub_start'] is not None:
        dates['sub_start'] = ee.Date(dates['sub_start'])
        try: dates['sub_start'].getInfo()
        except ee.ee_exception.EEException: raise ValueError(f"{Fore.RED}\033[1mInvalid sub period start date format; use YYYY-MM-dd.{Style.RESET_ALL}")
    dates['delta_sub'], dates['delta_period'] = NumChecker(dates['delta_sub'], 'Delta sub-period'), NumChecker(dates['delta_period'], 'Delta period')
    #Set to defaults if nothing entered
    if dates['delta_period'] is None: dates['delta_period'] = 1
    if dates['period'] is None: dates['period'] = 'year'
    #Check valid time units
    if dates['period'] not in time_units: raise ValueError(f'{Fore.RED}\033[1mPeriod must be a valid time unit.{Style.RESET_ALL}')

    #Global date checker
    DateChecker(sat, dates['start_date'], dates['end_date'])

    #Masking exception handling - some handled in the functions themselves
    masking['load_invert'], masking['save_raster'], masking['save_vector'] = FlagChecker(masking['load_invert'], 'mask loading inversion'), FlagChecker(masking['save_raster'], 'mask raster saving'), FlagChecker(masking['save_vector'], 'mask vector saving') 
    #Float conversions to allow for scientific notation read in 
    masking['upper_bound'], masking['lower_bound'] = NumChecker(masking['upper_bound'], "Mask upper bound"), NumChecker(masking['lower_bound'], "Mask lower bound")
    #Exception handling
    if masking['upper_bound'] is not None and masking['lower_bound'] is not None:
        if masking['upper_bound'] == 1e5 and masking['upper_bound'] == -1e-5: raise ValueError(f'{Fore.RED}\033[1mEither the lower or upper bound must be set to non-default values.{Style.RESET_ALL}')
        if masking['lower_bound'] >= masking['upper_bound']: raise ValueError(f"{Fore.RED}\033[1mLower bound of the mask must be below the upper bound.{Style.RESET_ALL}")
    if masking['load_overlap'] is None: masking['load_overlap'] = 1
    if masking['load_overlap'] != 1 and masking['load_overlap'] != 2: raise ValueError(f"{Fore.RED}\033[1mOverlap length in the mask loading must be set to 1 or 2.{Style.RESET_ALL}")
    if masking['load_name'] is not None and type(masking['load_name']) is not str: raise TypeError(f"{Fore.RED}\033[1mName of file to load a mask from must be a string{Style.RESET_ALL}")

    #Analyses exception handling
    #Index and quantitative analysis list handling 
    if type(quants['indices']) is str:
        inds = quants['indices'].strip().split()
        inds = [re.sub(',', '', ind) for ind in inds]
        quants['indices'] = inds
    if type(quants['quant_analyses']) is str:
        quans = quants['quant_analyses'].strip().split()
        quants['quant_analyses'] = [re.sub(',', '', quan) for quan in quans]
    #Area checks
    if quants['areas'] is not None:
        quant_areas = [ee.Geometry.Rectangle(coords) for coords in quants['areas'] if len(coords) > 2]
        quant_points = [ee.Geometry.Point(coords) for coords in quants['areas'] if len(coords) == 2]
        comb_areas = quant_areas+quant_points
        #Check that all these areas/points are valid
        for area in comb_areas:
            #Check that the areas are formatted correctly
            if area.getInfo()['type'] != 'Polygon' and area.getInfo()['type'] != 'Point': raise TypeError(f"{Fore.RED}\033[1mRegion for quantitative analysis not formatted as a valid area or point.{Style.RESET_ALL}")
            #Check that the area is bounded
            if area.isUnbounded().getInfo() == True: raise ValueError(f"{Fore.RED}\033[1mRegion for quantitative analysis is unbounded.{Style.RESET_ALL}")
            #Check that the areas sit within the AOI
            if aoi.contains(area).getInfo() == False: raise ValueError(f"{Fore.RED}\033[1mRegion for quantitative analysis is outside the overall AOI.{Style.RESET_ALL}") 
        
        quants.update({'areas': [quant_areas, quant_points]})
        quants['save_flag'] = FlagChecker(quants['save_flag'], 'Analysis area saving')
    # Bound checks
    quants['upper_bound'], quants['lower_bound'] = NumChecker(quants['upper_bound'], "Threshold analysis upper bound"), NumChecker(quants['lower_bound'], "Threshold analysis lower bound")

    # Type checks
    if quants['indices'] is None: quants['indices'] = []
    if quants['quant_analyses'] is None: quants['quant_analyses'] = []
    if type(quants['indices']) is not list: raise TypeError(f"{Fore.RED}\033[1mIndices list is not a list.{Style.RESET_ALL}")
    if type(quants['quant_analyses']) is not list: raise TypeError(f"{Fore.RED}\033[1mQuantitative analyses list is not a list.")
    if type(quants['bandName']) is not str and quants['bandName'] is not None: raise TypeError(f'{Fore.RED}\033[1mBand name for quantitative analyses must be a string.{Style.RESET_ALL}')
    #No analyses requested 
    if len(quants['indices']+quants['quant_analyses']) == 0: raise ValueError(f"{Fore.RED}\033[1mNo indices or analyses requested.{Style.RESET_ALL}")

    #Remove duplicate indices in the user lists - bit clunky 
    quants['indices'] = list(set(quants['indices']))
    #Checking that all requested analyses exist 
    bad_ans = [an for an in quants['indices'] if an not in Analysis_lib]
    if len(bad_ans) != 0: logger.warning(f"{Fore.YELLOW}Warning: Unrecognised analyses requested (and ignored): "+str(bad_ans)+". Valid analyses are "+str(Analysis_lib)+f'{Style.RESET_ALL}')
    quants['indices'] = [an for an in quants['indices'] if an not in bad_ans]

    quants['quant_analyses'] = list(set(quants['quant_analyses']))
    bad_ans = [an for an in quants['quant_analyses'] if an not in Quant_lib]
    if len(bad_ans) != 0: logger.warning(f"{Fore.YELLOW}Warning: Unrecognised analyses requested (and ignored): "+str(bad_ans)+". Valid analyses are "+str(Quant_lib)+f'{Style.RESET_ALL}')
    quants['quant_analyses'] = [an for an in quants['quant_analyses'] if an not in bad_ans]

    #Should compatibility handling be included here, in case of both making and defined vis?
    # Empty visualisation input handling - only defaults will be used in this case
    if vis is None: vis = {}
    # Checking for any custom visualisations and including these in the dictionary
    for ind in vis.keys():
        ind_vis = vis[ind]
        if ind_vis['make'] == 1:
            if ind_vis['floor'] is None: ind_vis.update({'floor': '#000000'})
            if ind_vis['ceiling'] is None: ind_vis.update({'ceiling': '#FFFFFF'})
            if ind_vis['n_steps'] is None: ind_vis.update({'n_steps': 10})
            cust_vis = VisMaker(ind_vis['vals'], ind_vis['colours'], ind_vis['steps_flag'], ind_vis['ramp_flag'], ind_vis['show'], names['PlotDir']+'/'+ind+'_Vis', ind_vis['floor'], ind_vis['ceiling'], ind_vis['n_steps'])
            vis.update({ind: cust_vis})
        else:
            ind_vis['vis_min'], ind_vis['vis_max'], ind_vis['palette'] = NumChecker(ind_vis['vis_min'], 'Visualisation minimum'), NumChecker(ind_vis['vis_max'], 'Visualisation maximum'), ColourListChecker(ind_vis['palette'])
            vis.update({ind: ind_vis})
            if ind_vis['show'] == 1:
                VisShower(ind_vis, names['PlotDir']+'/'+ind+'_Vis')

    return names, aoi, sat, process, dates, masking, quants, vis

