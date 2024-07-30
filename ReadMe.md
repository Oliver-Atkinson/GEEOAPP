# Cloud based processing flow for Google Earth Engine
Original work by Oliver Atkinson for Omanos Analytics (2022), repository:
https://gitlab.com/omanosanalytics/earth-engine-internship/

## Aims

Ideally, this project becomes something of a one stop shop for basic Google Earth Engine usage; loading in images according to the user specifications of region, satellite and time period, performing any requested analyses (including index calculations, masking, time period analyses such as area means) and directly saving any images and plots to Google Drive files.

A fair amount of this functionality is already present in the code, with a full explanation of the workings in the (as yet unwritten) manual. This ReadMe includes an outline of the package strcuture and a basic example of how to use the package to output some simple analysis of a given region. Lots of error checking goes on behind the scenes which is not outlined here.

## Setup

### Dependencies
* [Python 3](https://www.python.org/downloads/)
* [Google Earth Engine Python API](https://developers.google.com/earth-engine/guides/python_install) 

### Virtual Env

The requirements for the virtual enviornment used to run this project are given in GEEOAPP_CondaEnv.txt, which can be used to create a new virtual enviornment with ```conda create --name <env> --file  GEEOAPP_CondaEnv.txt```. NB: This was generated on an arm architecture mac, so things may require some tweaking. The list of packages is also likely to be overcomplete.

## Structure

A quick outline of what each of the pertinent python and yaml files contains and what their use is.

IndexFunctions.py - A library of index calculating functions. Indices with common structures, e.g. normalised difference indices (ND(Vegetation)I, ND(Water)I, ND(Snow)I) are collected together into single functions for ease. A library of RGB composites can also be called.

VisFunctions.py - Collects a number of functions for creating and processing visualisation parameters for image saving, including loading standard visualisations and capacity for user generation of new visualisations in a few different ways, with examples of these shown in VisFunctionality.pdf.

CollectionLoading.py - Holds the functions for requesting image collections from the Google Earth Engine database. The main loading function takes a satellite name, an area of interest and dates. Cloud masks can optionally applied, and collections filtered to only contain images with less than a given percentage of pixels cloud classified as cloud. If the requested satellite has no suitable images, other satellites that were operational at the same time will be checked and the user will be given the option to choose from alternatives (if any are available). In here are functions to generate image collections over a time series (e.g. every whole year, or every January of every year).

QuantFunctions.py - Stores the functionality for the quantitative analyses (mean, standard deviation, etc.) that can be performed on a given area or point, as well as the functions to plot time series graphs of the results of these analyses.

GEEOAPP.py - The central file, that imports all the functionality and information from all other files and is used to actually run the project. After imports, numerous error checks are performed on the user input from the run card. 
This code itself additionally contains key elements of the project, such as, all image masking capability (generating, saving and loading), and crucially, a function that brings all of the functionality together and serves to run the whole package.

GEE_Authenticator.py - A short piece of code that serves to check that the user is authenticated to use Google Earth Engine. 

GUI.py - Controls the graphical user interface, one of the possible input modes for the user; this code is only called if the user chooses not to supply a run card at run time.

Satellites.yml - Contains the library of currently supported satellites, storing the name of the Google Earth Engine database reference, the type of instrumentation (MSI, SAR, etc.), a dictionary of the bands in each image and what they physically correspond to (e.g for Sentinel 2, the band named "B3" corresponds to green), the resolution of the highest resolution band, the dates over which the satellite is/was operational, the cloud coverage property name, the cloud masking function name, and the scaling information.

Palettes.yml - Holds a large number of colour palettes, each of which is in the form of a list of hex strings, which can be used to create visualisations.

Visualisations.yml - Holds the standard visualisations for each index, with the form, storing the minimum and maximum values, and the name of the colour palette as it appears in Palettes.yml. 


## Usage

The package parameters are set in a run card by the user, and then read into the code at run time. This is an outline of how to perform a simple run of the code and save a true colour RGB image and a single band NDVI image. In this example we will look at Glasgow, creating median images of June through September for 2018-2022 from Sentinel 2. 

An example run card for this scenario is stored at `RunCards/Examples/Example1.yml`.

First, the identifying names for the run are set; which define where the resultant images will be saved and the prefix of their name. The area of interest is defined by the latitude and longitude of the four corners.

A satellite is chosen, which here is Sentinel 2; available options are held in the `Satellites.yml` file.

In the processing block, the collection filtering options are set; here the cloud mask is switched on and only images with less than 20% of cloudy pixels are passed. The `reducer` option tells the code how to collapse collections, of which there will be 5 here, into individual images. In this example, the median is used.

The dates define the period of interest, with the start and end dates setting the overall period. Here the collections are taken every year, between June and September, using the period and subperiod options.

No masking is performed here.

The only analyses of interest in this example are NDVI and RGB, which is the only input of import in the analyses block.

There is also no need to customise visualisations here.

Then run the code with a `python GEEOAPP.py Examples/Example1.yml` call in terminal.

The outputs to the terminal should be clear.

This will generate and save a total of 10 (5 RGB true colour, 5 NDVI) images, each with (hopefully) explanatory names indicating the dates of acquistion, nature of the image and satellite of origin. These will be saved directly to Google Drive in a folder named Example1, which will be created automatically if it does not already exist.

The package can also be run without providing a run card by calling `python GEEOAPP.py`. In this case, a graphical user interface will be launched and the user can input the various options here. These inputs are then checked for validity and a corresponding run card saved and automatically used if the inputs are indeed valid.