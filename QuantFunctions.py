#!/usr/bin/env python
"""A collection of functions for quantitative image analysis"""

# Imports
import ee
import numpy as np
import matplotlib.pyplot as plt

###############################################################
# To do:
#   Any more analyses to add
###############################################################

#To do image wide calculations, use e.g. im.select(band).reduceRegion(ee.Reducer.sum().unweighted(),aoi,res).get(band)

#Gives a mean in a given band over a given area 
def AreaMean(im, area, res, bandName):
    return im.reduceRegion(ee.Reducer.mean(),area,res).get(bandName).getInfo()

#Gives the standard deviation in a given band over a given area 
def AreaStdDev(im, area, res, bandName):
    return im.reduceRegion(ee.Reducer.stdDev(),area,res).get(bandName).getInfo()

#Gives the percentage of pixels in a given area that have values for a given band between the given upper and lower bounds
def AreaThreshold(im, area, res, bandName, upper_bound=1e5, lower_bound=-1e5):
    N_pix = ee.Number(im.reduceRegion(ee.Reducer.count(),area,res).get(bandName))
    lower = ee.Number(im.lte(upper_bound).reduceRegion(ee.Reducer.sum().unweighted(),area,res).get(bandName))
    upper = ee.Number(im.gte(lower_bound).reduceRegion(ee.Reducer.sum().unweighted(),area,res).get(bandName))
    n_g_pix = lower.add(upper).subtract(N_pix)
    g_perc = n_g_pix.divide(N_pix).multiply(100).getInfo()
    return g_perc

#Returns the value of a given band at a given point
def PointValue(im, point, res, bandName):
    return im.reduceRegion(ee.Reducer.first(),point,res).get(bandName).getInfo()

def TimePlotter(dates, y_vals, y_label, filename, name='', errs=np.zeros(1)):
    """Collects time series plotting capabilities, with flexibility for single and combined plots"""
        
    #Allows for handling of arrays with multiple data series
    n_data = 0
    if type(y_vals) is not list and y_vals.ndim != 1: n_data = y_vals.shape[0]

    #Plot parameters
    fig_width, fig_height, xlab_font, msize, linewidth = 12, 5, 16, 5, 1
    date_labels = [s_date[:4] for s_date in dates] #Currently just extracts year
    x_ticks = np.arange(len(dates))

    #If there are lots of dates, use only every other one on the axis and rescale plot parameters
    if len(date_labels) >= 10: 
        fig_width = len(date_labels)*np.sqrt(2)
        fig_height = fig_width/np.sqrt(5)
        xlab_font = fig_width
        msize *= len(date_labels)/10
        linewidth *= np.sqrt(3)
        date_labels = [date_labels[i] for i in np.arange(len(date_labels),step=2)]
        x_ticks = np.arange(len(date_labels)*2,step=2)

    #Doing the plotting
    plt.style.use('seaborn-whitegrid')
    plt.subplots(figsize=(fig_width,fig_height))

    #Dealing with the possibility of arrays and errors - somewhat clunky
    if n_data != 0:
        for i in range(n_data):
            plt.plot(y_vals[i],marker='o',label=name+' '+str(i+1),linewidth=linewidth, markersize=msize)
        plt.legend(fontsize=xlab_font)
        if errs.size == y_vals.size:
            x_points = np.arange(len(dates))
            for i in range(n_data):
                plt.fill_between(x_points, y_vals[i]-errs[i], y_vals[i]+errs[i],alpha=0.4)
    
    if n_data == 0: 
        plt.plot(y_vals,marker='o',linewidth=linewidth, markersize=msize)
        if len(errs) == len(y_vals):
            x_points = np.arange(len(dates))
            plt.fill_between(x_points, y_vals-errs, y_vals+errs,alpha=0.4)

    plt.xticks(ticks=x_ticks, labels=date_labels, fontsize=xlab_font, rotation=30)
    plt.yticks(fontsize=xlab_font)
    # plt.xlabel('Date', fontsize=xlab_font) #Maybe change to year or some other interval, possibly based on examination of the dates
    plt.ylabel(y_label, fontsize=xlab_font)
    plt.savefig(filename+'.pdf',bbox_inches='tight')

    return
