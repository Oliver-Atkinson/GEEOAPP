#!/usr/bin/env python
"""Collects the functions and, currently at least, the information required to load collections of images from GEE.
This includes the cloud masking and relevant scaling of images."""

#Imports
import ee
import datetime
import logging
import yaml
from colorama import Style, Fore, Back

logger = logging.getLogger(__name__)

# Reading in the satellite information
with open("Satellites.yml", "r") as stream:
    try: Satellite_Information = yaml.safe_load(stream)
    except yaml.YAMLError as exc: print(f'{Fore.RED}\033[1m'+str(exc))


def Sentinel2_CloudMask(im: ee.Image):
    """A cloud masking function for Sentinel 2 data"""
    qa = im.select('QA60') #QA60 band is the cloud mask

    #Bits 10 and 11 are clouds and cirrus, respectively
    cloudBitMask = 1 << 10
    cirrusBitMask = 1 << 11

    #Both flags should be set to zero, indicating clear conditions
    mask = qa.bitwiseAnd(cloudBitMask).eq(0) and (qa.bitwiseAnd(cirrusBitMask).eq(0))

    return im.updateMask(mask)

#---------Possible issue---------
#Onlys hold for the visible - different scaling and mask (QA_RADSAT.eq(0)) for thermal channels, which may not be present in all landsats
def Landsat_CloudMask(im: ee.Image):
    """A cloud masking function for Landsat data"""
    qa = im.select('QA_PIXEL') #QA_PIXEL band is the mask

    #Have the first 5 bits be clear (Fill, Dilated cloud, Cirrus, Cloud, Cloud shadow)
    cloud_mask = qa.bitwiseAnd(int('11111', 2)).eq(0)

    return im.updateMask(cloud_mask)

#For now, assume we want all possible bands, filtered for low cloud and a cloud mask applied
def CollectionGetter(sat: str, aoi: ee.Geometry, process_dict: dict, start_date: ee.Date, end_date: ee.Date):
    """Returns an image collection from a given satellite for a specified area and time period.
    
    Filtering for cloud and applying a cloud mask if requested, both of which vary by satellite, an image collection
    is generated between the start and end dates for the area of interest. If data matching the criteria cannot be found
    for the requested satellite, attempts to find data from other satellites of the same type that were operational
    during the requested time period, asking for user input on returning one of these instead if such data exists.
    
    Parameters:
    sat: str
        The name of the satellite to use data from
    aoi: ee.Geometry
        The area of interest 
    process_dict: dict
            A dictionary with the image processing parameters
                reducer: str - The name of the operation (currently median or mosaic) used to reduce collections to a single image
                cloud_cover: float - The value by which to filter for images with less than the given cloud coverage percentage
                cloud_mask: int [0,1] - Flag for application of the cloud mask
    start_date: ee.Date
        Start of the overall period
    end_date: ee.Date
        End of the overall period
        
    Returns:
    ee.ImageCollection
        The collection of images found - may be an empty collection, in which case a warning is raised
    str
        The satellite from which the images were found
    """
    
    #Actually does the business of getting a collection
    def ColMaker(sat):
        #Getting the relevant information from the satellite dictionary
        sat_inf = Satellite_Information[sat]
        col_name, col_type, col_cloud, col_scale_mult, col_scale_add, col_res = sat_inf['Dataset'], sat_inf['Type'], sat_inf['Cloud Property'], sat_inf['Scale Multiply'], sat_inf['Scale Add'], sat_inf['Resolution']

        #Scales the images in the collection, according to the requirments of the satellite dataset
        def Scaler(im):
            return im.multiply(col_scale_mult).add(col_scale_add)

        #Handling MSI and SAR seperately
        if col_type == 'MSI':
            #A bit clunky on the cloud masking, but straightforward options don't seem to work
            if process_dict['cloud_mask'] == 1: 
                col_mask = globals()[sat_inf['Mask']]   #Get the cloud masking function (no such requirement for SAR)

                #Making the image collection
                col = (ee.ImageCollection(col_name)     #Requesting the full dataset
                    .filterBounds(aoi)                  #Look at the area of interest
                    .filterDate(start_date, end_date)   #For the given date range
                    .filter(ee.Filter.lt(col_cloud, process_dict['cloud_cover'])) #Filter for images with less than the given cloud coverage
                    .map(col_mask)                      #Apply the cloud masking function if required
                    .map(Scaler))                       #Apply the scaling as defined above
            else:
                col = (ee.ImageCollection(col_name)     #Requesting the full dataset
                .filterBounds(aoi)                  #Look at the area of interest
                .filterDate(start_date, end_date)   #For the given date range
                .filter(ee.Filter.lt(col_cloud, process_dict['cloud_cover'])) #Filter for images with less than the given cloud coverage
                .map(Scaler))                       #Apply the scaling as defined above                            

        #Plenty of other options can be done here - ascending or descending pass, instrumentMode and resolution options
        if col_type == 'SAR':
            col = (ee.ImageCollection(col_name)     #Requesting the full dataset
                .filterBounds(aoi)                  #Look at the area of interest
                .filterDate(start_date, end_date)   #For the given date range
                .map(Scaler))                       #Apply the scaling as defined above

        #Stamping relevant information onto the collection 
        col = col.set({'Satellite': sat,
                'StartDate': start_date.format('YYYY-MM-dd'),
                'EndDate': end_date.format('YYYY-MM-dd'),
                'Resolution': col_res})
    
        return col
    
    col = ColMaker(sat)

    #If the collection is empty, send a warning about this and check for possible alternatives
    if col.size().getInfo() == 0:
        logger.warning(f"{Fore.YELLOW}Warning: There is no "+sat+" data for "+str(start_date.format('YYYY-MM-dd').getInfo())+" - "+str(end_date.format('YYYY-MM-dd').getInfo())+f" for the region selected that satisfies the requested criteria.{Style.RESET_ALL}")
        print("Attempting to find datasets which do have images for the requested time and place...")
                
        #First, get a list of satellites of the same type that were operational during the requested period
        #List of satellites of the same type (remove the original satellite to avoid repeats)
        sat_type = Satellite_Information[sat]['Type']
        good_type_sats = [gsat for gsat in Satellite_Information.keys() if Satellite_Information[gsat]['Type'] == sat_type and gsat != sat]
        #For all the same type satellites, check if the satellite was operational during the requested dates
        bad_date_sats = []
        for tsat in good_type_sats:
            #Checking that the satellite was operational during the dates requested - stolen from DateChecker, quicker than running full DateChecker
            sat_dates = Satellite_Information[tsat]['Operative']
            sat_start, sat_end = sat_dates[0], sat_dates[1]
            if sat_dates[1] == 'Present': sat_end = ee.Date(str(datetime.date.today())) #If the satellite is still active, set the end date to today - not always accurate as data not up to the exact date
            sat_range = ee.DateRange(ee.Date(sat_start), ee.Date(sat_end))
            #Check if the dates are valid
            req_range = ee.DateRange(start_date, end_date)
            if req_range.intersection(sat_range).isEmpty().getInfo() == True: bad_date_sats.append(tsat)

        good_sats = [gsat for gsat in good_type_sats if gsat not in bad_date_sats]
        
        #Dictionary to hold possible alternative collections, allowing for choice by the user
        alt_cols = {}

        #For every satellite in the appropriate list, try to find a collection for the given area/times
        for test_sat in good_sats: 
            test_col = ColMaker(test_sat) 
            N_ims = test_col.size().getInfo()
            if N_ims != 0:
                alt_cols.update({test_sat: test_col})
                #If a collection has any images in it, print a message saying so, with the resolution of the collection
                if N_ims == 1: print("1 suitable image found in the "+test_sat+" dataset ("+str(Satellite_Information[test_sat]['Resolution'])+"m resolution).")
                else: print(str(N_ims)+" suitable images found in the "+test_sat+" dataset ("+str(Satellite_Information[test_sat]['Resolution'])+"m resolution).")
        
        #List of full satellite names that do have data
        alt_sats = list(alt_cols.keys())

        #If no images found anywhere, raise a warning (maybe upgrade to an error?), and return the original (empty) collection
        if len(alt_sats) == 0: 
            logger.warning(f"{Fore.YELLOW}Warning: No images found in any dataset - an empty collection will be returned, likely leading to errors...{Style.RESET_ALL}")
            return col, sat
        
        #If only one other satellite has available data, ask a straight up or down
        if len(alt_sats) == 1:
            y_n = input("This is the only dataset with available data; use this instead? (y/n)\n")
            if y_n == "y" or y_n == "yes" or y_n == "Yes": return alt_cols[alt_sats[0]], alt_sats[0]
            else: 
                logger.warning(f"{Fore.YELLOW}Warning: An empty collection will be returned, likely leading to errors...{Style.RESET_ALL}")
                return col, alt_sats[0]

        #List of satellite names for easier user input, e.g. Sentinel 2 -> S2 
        name_list = [sat_name[0]+sat_name[-1] for sat_name in alt_sats]
        name_dict = dict(zip(name_list, alt_sats)) #Dictionary to map abbreviations to full names
        #Create a string to print detailing the user's options
        p_string = "Enter "
        for n in range(len(name_list)):
            p_string += name_list[n]+" for "+alt_sats[n]+", "
        p_string += "or any other input for the empty collection (not advised):\n"

        #Allowing the user to choose an alternative satellite which does have some data for the aoi at the chosen time
        print("Choose an alternate satellite dataset to use, or proceed with an empty collection (not advised).")
        alt_col_choice = input(p_string)
        if alt_col_choice not in name_list: 
            logger.warning(f"{Fore.YELLOW}Warning: An empty collection will be returned, likely leading to errors...{Style.RESET_ALL}")
            return col
        col = alt_cols[name_dict[alt_col_choice]]
        return col, name_dict[alt_col_choice]

    return col, sat


def TimeSeries(aoi: ee.Geometry, sat: str, process_dict: dict, start_date: ee.Date, end_date: ee.Date, delta_interval: float, interval: str):
    """Creates a collection of images, each from an interval within a given time period, e.g every 3 months
    
    By first creating a series of dates, a collection of images is generated, where each image is a [median] over
    a given interval, such as every month, within the overall start and end dates.
    
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
    start_date: ee.Date
        Start of the overall period
    end_date: ee.Date
        End of the overall period
    delta_interval: float
        Number of interval time units to step, e.g. every delta_interval days/months/years  
    interval: str
        The time unit to step, e.g. every n interval

    Returns:
    list
        A list of images, each of which is a over the interval specified, within the start and end dates given
    """

    print('Finding images from each '+'{0:.3g}'.format(delta_interval)+' '+interval+' period between '+start_date.format('YYYY-MM-dd').getInfo()+' and '+end_date.format('YYYY-MM-dd').getInfo()+'...')

    #Find the number of individual date ranges to examine
    window = end_date.difference(start_date, interval)
    n_intervals = window.divide(delta_interval).int().getInfo()
    #If there's only one image, which is the minimum, raise a warning
    if n_intervals == 0: logger.warning(f"{Fore.YELLOW}Warning: Dates are set such that only a single collection will be returned.{Style.RESET_ALL}")
    
    #Creating a list of dates at given intervals within the overall window, with EE functionality
    #Function to iterate over
    def DateStepper(current_elem, previous_list):
        prev_date = ee.Date(ee.List(previous_list).get(-1))      #Get the last value of the list
        next_date = prev_date.advance(delta_interval, interval)  #Advance it delta_interval intervals
        return ee.List(previous_list).add(next_date)             #Return a list with the new date added
   
    nDates = ee.List.sequence(1, n_intervals) #Sequential list to run the iteration on (analogue of range(n_intervals))
    dates = ee.List(nDates.iterate(DateStepper, [start_date])).add(end_date) #Do the iteration, making a list of start_date, and then append the end_date
   
    #Warning about final period being curtailed
    fin_forward = ee.Date(dates.get(-2)).advance(delta_interval, interval)
    if end_date.difference(fin_forward, 'day').getInfo() < 0:
        logger.warning(f'{Fore.YELLOW}Warning: The final time period requested has an end date ('+fin_forward.format('YYYY-MM-dd').getInfo()+f') beyond the overall end date, and so has been curtailed to end at the overall end date.{Style.RESET_ALL}')

    #For each interval, get the collection, and then add to a list
    col_list = [CollectionGetter(sat, aoi, process_dict, ee.Date(dates.get(n)), ee.Date(dates.get(n+1))) for n in range(n_intervals+1)]

    #Return a list of image collections
    return col_list

#This code will by nature be similar and repetitive to the TimeSeries function in parts at least
def SubperiodSeries(aoi: ee.Geometry, sat: str, process_dict: dict, start_date: ee.Date, end_date: ee.Date, delta_sub: float, subperiod: str, sub_start: ee.Date = ee.Date(0), delta_period: float = 1, period: str = 'year'):
    """Creates a collection of images for given sub-periods within the overall period, e.g. every January from 2010-2020
    
    By first creating a series of non-overlapping dates, iterating over this to generate collections of images, projected
    down to one image per sub-period by taking a [median], a collection of images from the entire start to end date is generated.
    
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
    start_date: ee.Date
        Start of the overall period
    end_date: ee.Date
        End of the overall period
    delta_sub: float
        Length of sub-period time to sample, in units of subperiod, e.g. a duration of delta_sub days/months/years  
    subperiod: str
        Time unit for the sub-period, e.g. delta_sub subperiods
    sub_start: ee.Date (Optional)
        Start of the sub-period to sample (must be within the overall start and end dates), defaults to the overall start_date
    delta_period: float (Optional)
        Number of time units to step between sub-periods, e.g. sample every delta_period days/months/years, defaults to 1  
    period: str (Optional)
        The time unit to step between sub-periods, e.g. every n period, defaults to year

    Returns:
    list
        A list of images, each of which is over a single sub-period, within the start and end dates given
    """

    #If no sub_start is passed, assume it's the same as the overall start date
    if sub_start == ee.Date(0): sub_start = start_date
    
    #Exception handling
    #Check the sub period is less than the iterative period - bit ungainly
    if ee.DateRange(start_date.advance(delta_period, period), start_date.advance(delta_sub, subperiod)).isEmpty().getInfo() == False: 
        logger.warning(f"{Fore.YELLOW}Warning: Sub-period longer than the period - this is fine for a moving average but is not generally ideal.{Style.RESET_ALL}")
    

    #Check sub period start lies within the larger period
    req_range = ee.DateRange(start_date, end_date)
    if req_range.contains(sub_start).getInfo() == False: raise ValueError(f"{Fore.RED}\033[1mSub period start date not within the overall period.")
    
    print('Finding images from each '+'{0:.3g}'.format(delta_sub)+' '+subperiod+' period starting from '+sub_start.format('YYYY-MM-dd').getInfo()+' at '+'{0:.3g}'.format(delta_period)+' '+period+' intervals between '+start_date.format('YYYY-MM-dd').getInfo()+' and '+end_date.format('YYYY-MM-dd').getInfo()+'...')

    #Find the number of individual date ranges to examine
    n_periods = end_date.difference(start_date, period).getInfo()
    n_intervals = int(n_periods/delta_period) + 1 #The +1 accounts for the int call rounding
    #If there's only one image, which is the minimum, raise a warning
    if n_intervals == 1: logger.warning(f"{Fore.YELLOW}Warning: Dates are set such that only a single collection will be returned.{Style.RESET_ALL}")

    #Creating a list of pairs of dates, one pair (start and end) for each subperiod
    dates = []
    for n in range(n_intervals):
        sub_0 = sub_start.advance(n*delta_period, period) #Move the start date one full period each time
        sub_1 = sub_0.advance(delta_sub, subperiod)       #Advance the end date one subperiod from the start of the subperiod
        dates.append([sub_0, sub_1])

    #Check the sub period dates are valid, i.e. sub period end dates are before the overall end date
    #Using ee.Filter.dateRangeContains might be more natural, but this is complicated
    end_dates = ee.List([date[1].millis() for date in dates])             #Get an EE list of end dates as numbers
    filt_end = end_dates.filter(ee.Filter.lte('item', end_date.millis())) #Filter the list to give end_dates before the overall end date
    last_good_ind = end_dates.indexOf(filt_end.get(-1)).add(1).getInfo()  #We only care about the index of the last end_date (+1 for slicing) before the overall end date 
    good_dates = dates[:last_good_ind]                                    #Slicing the dates list to remove pairs with bad end dates 

    #Include warning about removed dates - could specify the actual date periods removed if this is worth doing?
    if len(good_dates) != len(dates): logger.warning(f'{Fore.YELLOW}Warning: The last '+str(len(dates)-len(good_dates))+f' sub-period(s) had end dates beyond the overall end date, and so have been ignored.{Style.RESET_ALL}')

    # For each interval, get the collection and add it to a list
    col_list = [CollectionGetter(sat, aoi, process_dict, date_pair[0], date_pair[1]) for date_pair in good_dates]

    #Return the list of collections
    return col_list
