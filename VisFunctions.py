#!/usr/bin/env python

"""Collects a number of functions for creating and processing visualisation parameters for image saving."""

import logging
import yaml
from colorama import Style, Fore, Back
import numpy as np
import matplotlib.pyplot as plt
import re

from colour import Color

logger = logging.getLogger(__name__)

###############################################################
# To do:
#   Expand visualisation library
#   Error handling - any other fail cases?
#   May be repetitive/redundant calls, particularly on error handling and Maker/Loader 
###############################################################


# Reading in palettes
with open("Palettes.yml", "r") as stream:
    try: Palettes = yaml.safe_load(stream)
    except yaml.YAMLError as exc: print(f'{Fore.RED}\033[1m'+str(exc))

# Reading in visualisations
with open("Visualisations.yml", "r") as stream:
    try: Visualisations = yaml.safe_load(stream)
    except yaml.YAMLError as exc: print(f'{Fore.RED}\033[1m'+str(exc))

#################################

def ClampRGB(n, lower:int=0, upper:int=255): 
    """Clamps a number n to be between lower and upper bounds"""
    return max(lower, min(n, upper))

def HexCheck(hex: str):
    """Checks if a string is a valid hex code"""
    #Must be a string
    if type(hex) is not str:
        return False

    #Starts with '#', 3 or 6 other characters, each either 0-9, a-z, A-Z
    match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hex)
    if match: return True
    else: return False

def hex_to_rgb(value: str):
    """Return (red, green, blue) array for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return np.array([int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3)])

def rgb_to_hex(rgb):
    """Return a color as #rrggbb from a given list or array of RGB values, either as floats 0-1 or ints 0-255."""
    #Exception handling
    #Type check
    if type(rgb) is not type(np.zeros(1)) and type(rgb) is not list: raise TypeError(f"{Fore.RED}\033[1mRGB colours must be formatted as a list or array")
    #Must be of length 3
    if len(rgb) != 3: raise ValueError(f"{Fore.RED}\033[1mRGB format must be a list or array of length 3")

    #-------Possible issue-------
    #Extremal cases of 0 and 1 in RGB floats will be handled badly - i.e. assumed to be valid ints within 0-255
    #If 0 < rgb values < 1, scale by 255 - i.e. 0-255 as the default format, and allow for conversion otherwise
    #Clamp any values outwith 0-255 and force conversion to ints
    rgb = [int(val*255) if 0 < val < 1 else int(ClampRGB(val)) for val in rgb]

    return '#%02x%02x%02x' % (rgb[0], rgb[1], rgb[2])

def ColourListChecker(palette):
    """Checks if a list of colours is valid, converting RGB to hex if required and returning a list of hex strings"""
    #Type check
    if type(palette) is not str and type(palette) is not list: raise TypeError(f"{Fore.RED}\033[1mColour list to make a palette from must be a palette name (string) or a list of colours (hex strings or [R,G,B] values).")

    #If a string is passed, look up the palette in Palettes.yml and do no further checks
    if type(palette) is str: 
        try: palette = Palettes[palette]
        except KeyError: raise ValueError(f"{Fore.RED}\033[1mNo colour list called "+palette+" exists - see Palettes.yml for valid palettes.")
        return palette
        
    #Palette must be a list at this point, now check the entries
    good_palette = [] #New list to allow for conversions and changes as a result of the error handling
    for col in palette:
        #Add '#' if the colour is a non-empty string that doesn't already start with '#' already (some in the Palettes.yml don't)
        if type(col) is str and len(col) != 0 and col[0] != '#': col = '#'+col
        #If not a string, assume the colour is in an RGB format and convert to hex
        if type(col) is not str:
            col = rgb_to_hex(col) 

        #Check if the colour is a valid hex
        if HexCheck(col) == False: raise ValueError(f"{Fore.RED}\033[1mInvalid hex string '"+col+"' in the colour list.")

        #Append the resulting checked colour to the valid palette
        good_palette.append(col)

    return good_palette

def ColourRamp(colours: list, n_steps: int = 10):
    """Creates a ramp of colours between the colours given in the colours parameter;
    (colours[0] -> colours[1] -> ... -> colours[-1]) with a total of n_steps colours in the final ramp list.
    
    Parameters:
    colours: list
        A list of colours (each hex or [R,G,B]), with each colour being a start point of a ramp
    n_steps: int (Optional)
        The total number of steps in which to go from the first colour to the last. Defaults to 10.
    
    Returns:
    Ramp: list
        A list of hex strings, of length steps, that ramps between each colour pair in colours
    """
    #Exception handling
    if type(n_steps) is not int: raise TypeError(f"{Fore.RED}\033[1mNumber of steps for the colour ramp must be an integer.")
    if len(colours) > n_steps: raise ValueError(f"{Fore.RED}\033[1mMore colours than steps requested - would result in missing colours.")
    colours = ColourListChecker(colours)

    #Changing the given list to a list of colour objects
    Cols = [Color(col) for col in colours]

    #Number of ramps required is one less than the number of colours
    N_ramps = len(Cols)-1 

    #The number of steps per ramp required is slightly more involved to ensure the final ramp has len == steps
    ramp_step_arr = np.ones(N_ramps)*int(n_steps/N_ramps) 
    #int(steps/N_ramps) alone will lead to an undershoot, due to rounding errors and the removal of ramp end/start pairs when joining multiple ramps
    undershoot = (n_steps % N_ramps) + N_ramps-1
    #Distribute the undershoot amongst the ramp step array until it's taken care of
    while undershoot >= N_ramps:
        ramp_step_arr += 1
        undershoot -= N_ramps
    #Slightly stretch the first ramp(s) as required
    for i in range(undershoot):
        ramp_step_arr[i] += 1

    # Sanity check on number of steps per ramp
    if min(ramp_step_arr) <= 3: logger.warning(f"{Fore.YELLOW}Warning: Low number of steps per ramp colour ("+str(int(min(ramp_step_arr)))+f") - may result in a non-smooth ramp or missing colours.{Style.RESET_ALL}")
    
    #Creating the ramp, first as a list of individual ramps, then concatenating
    Ramp_lists = [list(Cols[n].range_to(Cols[n+1], int(ramp_step_arr[n]))) for n in range(N_ramps)] #List of colour ramps
    Ramp = []
    for i in range(N_ramps-1):
        Ramp += Ramp_lists[i][:-1] #The last colour in each ramp is the start of the next - to avoid repeats we don't include it
    Ramp += Ramp_lists[-1] #The full final ramp should be included (hence the -1 in range(N_ramps-1)), to ensure the final colour ends the overall ramp
    Ramp = [col.hex for col in Ramp] #Converting to hex strings

    return Ramp

def VisFromSteps(steps: list, colours: list, floor='#000000', ceiling='#FFFFFF', ramp_flag=0):
    """Create a visualisation from a series of steps, with a palette assigned to intervals i.e. if between steps[i] and steps[i+1], be colours[i]
    
    Parameters
    steps: list
        A list of floats corresponding to the end points of intervals. 
    colours: list
        The colour list used to create the palette.  This must have length one less than the length of steps. The ith
        interval in steps (steps[i] to steps[i+1]) will correspond to colours[i] in the visualisation.
    floor: str (Optional)
        The colour that will be used for values below steps[0] - defaults to black.
    ceiling: str (Optional)
        The colour that will be used for values above steps[-1] - defaults to white. 
    ramp_flag: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on generates a colour ramp between each colour in the colours parameter and uses this as the palette.
        
    Returns
    vis: dict
        A visualisation dictionary {vis_min: <float>, vis_max: <float>, palette: list of hex strings}
    """
    #Exception handling
    #Type handling
    if type(steps) is not list: raise TypeError(f"{Fore.RED}\033[1mThe steps to be used in creating a visualisation from steps must be passed as a list of floats.")
    for val in steps:
        try: float(val)
        except: raise TypeError(f"{Fore.RED}\033[1mInvalid number '"+str(val)+"' in the steps visualisation list.")
    # Steps should be an ascending list 
    if sorted(steps) != steps: raise ValueError(f"{Fore.RED}\033[1mThe steps visualisation list must be in ascending order.")
    
    #If a colour ramp is requested, using this instead
    if ramp_flag == 1:
        colours = ColourRamp(colours, len(steps)-1) #-1 as we want the number of intervals
    
    # Number of colours should be one less than the number of steps, i.e. equal to the number of intervals
    if len(colours)+1 != len(steps): raise ValueError(f"{Fore.RED}\033[1m""There must be one fewer colour than number of steps when creating a visualisation from steps; \
there were """+str(len(steps))+""" steps and """+str(len(colours))+""" colours passed.""")
    #Checking the colour list
    colours = ColourListChecker(colours)

    #Returns the colour from the palette that corresponds to the value passed
    def Val2Col(it_val):
        for j in range(len(steps)-1):
            if it_val <= steps[j]:
                return colours[j]            
            elif steps[j] < it_val <= steps[j+1]:
                return colours[j+1]
   
    #Minimum step for fine-graining the palette - find min difference, but then adjust down to match required step (for e.g. 0,0.1,0.25)
    step_0 = min([steps[i+1]-steps[i] for i in range(len(steps)-1)])
    #If the maximum remainder from dividing a step by the minimum step distance is less than the min step seperation, redefine the min_step
    if max([round((step%step_0)/step_0,5) for step in steps]) < 1:
        min_step = min([steps[i+1]-steps[i]-step_0 for i in range(len(steps)-1) if steps[i+1]-steps[i]-step_0 > 1e-10])
    else: min_step = step_0

    #Iterate over the step range, choosing a colour to add for each step
    colours = [floor]+colours #For values below the minimum, the floor colour will be used in the visualisation
    palette = []
    for i in np.arange(min(steps), max(steps)+1e-5, min_step):
        palette.append(Val2Col(round(i,5))) #Rounding to avoid slight overstepping
    palette.append(ceiling) #For values above the maximum, the ceiling colour will be used in the visualisation

    vis = {'vis_min': steps[0], 'vis_max': steps[-1], 'palette': palette}
    
    return vis

def VisShower(vis: dict, savename: str = "VisualisationFig", steps: list = []):
    """Produces a plot to show the visualisation scale.
    
    Parameters:
    vis: dict
        A visualisation dictionary {vis_min: <float>, vis_max: <float>, palette: list of hex strings}
    savename: str
        The filename by which to save the plot of the visualisation, defaults to VisualisationFig.pdf
    steps: list
        A list of floats, of possibly uneven intervals, used when showing a visualisation produced by the VisFromSteps function
     """
    
    vis_min, vis_max = vis['vis_min'], vis['vis_max']
    #Flexible palette type handling
    if type(vis['palette']) is str: palette = Palettes[vis['palette']]
    elif type(vis['palette']) is list: palette = vis['palette']

    #Checking if min and max for the visualisation are valid
    if vis_min >= vis_max: raise ValueError(f"{Fore.RED}\033[1mMinimum value greater than the maximum for the visualisation (vis_min >= vis_max)")
    #Check the palette is valid
    palette = ColourListChecker(palette)

    #Some setup values
    #If creating a plot from a series of (possibly uneven) steps, things should be handled differently
    if len(steps) != 0: 
        #Applying the same process to find the appropriate min step as in the VisFromSteps method
        step_0 = min([steps[i+1]-steps[i] for i in range(len(steps)-1)])
        if max([round((step%step_0)/step_0,5) for step in steps]) < 1:
            step = min([steps[i+1]-steps[i]-step_0 for i in range(len(steps)-1) if steps[i+1]-steps[i]-step_0 > 1e-10])
        else: step = step_0

        #The labels and ticks follow straight from the steps themselves in this case
        x_labels = ["0" if -1e-5 < val < 1e-5 else "{0:.3g}".format(val) for val in steps]
        x_ticks = steps    
    else: 
        step = (vis_max-vis_min)/(len(palette)-2)
        x_labels = ["0" if -1e-5 < val < 1e-5 else "{0:.3g}".format(val) for val in np.arange(vis_min, vis_max+step,step)]    
        x_ticks = np.arange(vis_min, vis_max+step,step)
    rng = np.arange(vis_min-step, vis_max+step+1e-5,step)

    #Some formatting options to ensure the figure displays well (up to ~20 steps)
    fig_width, fig_height, xlab_font = 12, 5, 16
    if 11 < len(x_labels) <= 21:
        fig_width = len(x_labels)*np.sqrt(2)
        fig_height = fig_width/np.sqrt(5)
        xlab_font = fig_width
    #If there are too many steps, don't show the x_labels and use default sizing
    if len(x_labels) > 21:
        x_labels = []

    #Creating the plot
    fig, ax = plt.subplots(figsize=(fig_width,fig_height))
    for i in range(len(rng)-1):
        ax.axvspan(rng[i], rng[i+1], color=palette[i])
    plt.xticks(ticks=x_ticks, labels=x_labels, fontsize=xlab_font)
    plt.yticks(ticks=[])
    plt.xlim([vis_min-step,vis_max+step])
    plt.savefig(savename+'.pdf',bbox_inches='tight')

    return

def VisMaker(vals: list, colours: list, steps_flag: int=0, ramp_flag: int=0, show: int=0, show_name: str='VisualisationFig', floor: str='#000000', ceiling: str='#FFFFFF', n_steps: int=10):
    """Creates a visualisation dictionary, using a series of steps, and/or a colour ramp
    created from a given list of colours, with settings based on user input.

    Parameters
    vals: list
        The values for the visualisation, as a list of steps in the visualisation from steps case, or list [min, max] in the ramp case
    colours: list
        The colours from which to form a palette, one for each interval in the case of only visualisation from steps, or the colours
        from which to create a ramp if the ramp_flag is on
    steps_flag: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on generates the visualisation using the VisFromSteps method. Defaults to 0 (off)
    ramp_flag: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on generates a colour ramp, which is used as the palette. Defaults to 0 (off)
    show: int [0,1] (Optional)
        Flag (0/1 for off/on) - if on saves a plot of the visualisation
    show_name: str (Optional)
        Filename for the visualisation display figure
    floor: str (Optional)
        The colour that will be used for values below vals[0] - defaults to black.
    ceiling: str (Optional)
        The colour that will be used for values above vals[-1] - defaults to white. 
    n_steps: int (Optional)
        The total number of steps in which to go from the first colour to the last, if the visualisation is created from only a colour ramp. 
        Defaults to 10.

    Returns
    vis: dict
        A visualisation dictionary {vis_min: <float>, vis_max: <float>, palette: list of hex strings}
    """
    
    #Exception handling
    #Flag value checks
    if steps_flag == 0 and ramp_flag == 0 and show == 0: raise ValueError(f"{Fore.RED}\033[1mNo flags switched on in visualisation maker - at least one of steps_flag, ramp_flag and show must be set to 1.")
    if steps_flag != 0 and steps_flag != 1: raise ValueError(f"{Fore.RED}\033[1mInvalid steps flag value in visualisation maker - only 0 and 1 accepted.")
    if ramp_flag != 0 and ramp_flag != 1: raise ValueError(f"{Fore.RED}\033[1mInvalid ramp flag value in visualisation maker - only 0 and 1 accepted.")
    if show != 0 and show != 1: raise ValueError(f"{Fore.RED}\033[1mInvalid show flag value in visualisation maker - only 0 and 1 accepted.")
    #Floor and ceiling hex checks
    if HexCheck(floor) != True: raise ValueError(f"{Fore.RED}\033[1mInvalid hex string for palette floor: "+str(floor))
    if HexCheck(ceiling) != True: raise ValueError(f"{Fore.RED}\033[1mInvalid hex string for palette ceiling: "+str(ceiling))
    #Checks on the colour list are handled in each of the individual functions - somewhat repetitively in some cases

    #If requesting a visualisation from steps, run this, passing the ramp_flag as an argument
    if steps_flag == 1: 
        vis = VisFromSteps(vals, colours, floor=floor, ceiling=ceiling, ramp_flag=ramp_flag)
        #If requested, save an image showing the visualisation scale, taking into account the steps
        if show == 1: VisShower(vis, show_name, vals)

    #If a colour ramp only is requested, generate this
    elif ramp_flag == 1:
        #Get the palette from the colour ramp, adding in the floor and ceiling values 
        palette = ColourRamp(colours, n_steps)
        palette = [floor]+palette
        palette.append(ceiling)
        
        #When creating a visualisation from a colour ramp only, the val list should only be the min and max value ([min,max])
        if type(vals) is not list: raise TypeError(f"{Fore.RED}\033[1mThe values used to create a visualisation with a colour ramp must be passed as a list of floats.")
        for val in vals:
            try: float(val)
            except: raise TypeError(f"{Fore.RED}\033[1mInvalid number '"+str(val)+"' in the colour ramp values list.")
        if len(vals) <= 1: raise ValueError(f"{Fore.RED}\033[1mThe list of values used to create a visualisation with a colour ramp should be just the min and max values, i.e. length 2, not "+str(len(vals)))
        if len(vals) > 2: logger.warning(f"{Fore.YELLOW}Warning, the list of values used to create a visualisation with a colour ramp should be just the min and max values. The first and last values in "+str(vals)+f" will be used.{Style.RESET_ALL}")

        vis = {'vis_min': vals[0], 'vis_max': vals[-1], 'palette': palette}
        #If the flag for showing the visualisation is on, save an image showing the visualisation scale
        if show == 1: VisShower(vis, show_name)

    #If the only request is to show the visualisation, do this
    else:
        vis = {'vis_min': vals[0], 'vis_max': vals[-1], 'palette': ColourListChecker(colours)}
        if show == 1: VisShower(vis, show_name)

    return vis

def VisLoader(vis, index: str):
    """Returns a valid visualisation dictionary for a given index, taking default settings if no input is given, otherwise
    looking up a palette or visualisation as requested."""

    #Looking up a visualisation if passed a visualisation name
    if type(vis) is str:
        try: vis = Visualisations[vis]
        except KeyError: raise ValueError(f"{Fore.RED}\033[1mNo visualisation called "+vis+" exists - see Visualisations.yml for valid visualisations.")
        
    #Weird stuff going on with definitions being overwritten - a bodge to correct this
    if type(vis['palette']) is str: vis = {'vis_min': vis['vis_min'], 'vis_max': vis['vis_max'], 'palette': Palettes[vis['palette']]}

    # If vis dictionary is empty, return the default visualisation for the index, looking up the palette in the process
    if len(list(vis.keys())) == 0:
        vis = Visualisations[index+'_default']
        vis['palette'] = Palettes[vis['palette']]
        return vis

    #Values from the visualisation dictionary, split out for checking
    vis_min, vis_max, palette = vis['vis_min'], vis['vis_max'], vis['palette']

    # Checks on types for non-zero input
    if type(vis_min) is not float and type(vis_min) is not int and type(vis_min) is not list: raise TypeError(f"{Fore.RED}\033[1mMinimum value to be used in the visualisation must be a float or int.")
    if type(vis_max) is not float and type(vis_max) is not int and type(vis_max) is not list: raise TypeError(f"{Fore.RED}\033[1mMaximum value to be used in the visualisation must be a float or int.")
    #Min and max for the visualisation are invalid
    if vis_min >= vis_max: raise ValueError(f"{Fore.RED}\033[1mMinimum value greater than the maximum for the visualisation (vis_min >= vis_max)")
    #Palette checking
    palette = ColourListChecker(palette)

    return vis
