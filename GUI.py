"""A simple GUI for GEEOAPP input data entry"""

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from RunCardReader import *

#Global variables
NDIs_lib = ['NDVI','NDWI','NDMI','NDSI','NDRE','GNDVI', 'NBR']#,'NBR' - no NBR vis yet, but it can be calculated 
Ratm1s_lib = ['RECI','GCI']
SAVIs_lib = ['SAVI','OSAVI','MSAVI']
Misc_lib = ['SIPI','EVI','VARI','ARVI', 'ARI', 'MARI', 'PSRI', 'OSI', 'BSI']
Index_lib =  NDIs_lib+Ratm1s_lib+SAVIs_lib+Misc_lib
RGB_lib = ['RGB', 'RedVegetation1','RedVegetation2', 'SWIRComposite1', 'SWIRComposite2', 'Agriculture','Bathymetric1','Bathymetric2','Bathymetric3','Bathymetric4','Bathymetric5','Bathymetric6','Urban']
RGBCalc_lib = ['Neon']
Analysis_lib = Index_lib+RGB_lib+RGBCalc_lib

Quant_lib = ['AreaMean', 'AreaStdDev', 'AreaThreshold', 'PointValue']



#Initialise root window
root = tk.Tk()
root.title("GEEOAPP")
tabControl = ttk.Notebook(root)

#Setting up the tabs
main_tab = ttk.Frame(tabControl)
mask_tab = ttk.Frame(tabControl)
analysis_tab = ttk.Frame(tabControl)

#Tab names
tabControl.add(main_tab, text ='Main') #Better name here...
tabControl.add(mask_tab, text ='Masking')
tabControl.add(analysis_tab, text="Analysis")

#Expanding
tabControl.pack(expand = 1, fill ="both")

#----------------------------------------
#               Main tab
#----------------------------------------

# Names
name_info_frame = tk.LabelFrame(main_tab, text="File naming information")
name_info_frame.grid(row= 0, column=0, padx=20, pady=10)

run_name_label = ttk.Label(name_info_frame, text="Run Name").grid(row=0, column=0)
imsave_label = ttk.Label(name_info_frame, text="Image saving directory").grid(row=0, column=1)
masksave_label = ttk.Label(name_info_frame, text="Mask saving directory").grid(row=0, column=2)
plotsave_label = ttk.Label(name_info_frame, text="Plot saving directory").grid(row=0, column=3)

run_name_entry = ttk.Entry(name_info_frame)
run_name_entry.grid(row=1, column=0)
imsave_name_entry = ttk.Entry(name_info_frame)
imsave_name_entry.grid(row=1, column=1)
masksave_name_entry = ttk.Entry(name_info_frame)
masksave_name_entry.grid(row=1, column=2)
plotsave_name_entry = ttk.Entry(name_info_frame)
plotsave_name_entry.grid(row=1, column=3)


for widget in name_info_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)
    
# AOI
aoi_frame = tk.LabelFrame(main_tab, text="Area of interest")
aoi_frame.grid(row=1, column=0, sticky="news", padx=20, pady=10)

min_long_label = tk.Label(aoi_frame, text="Minimum Longitude").grid(row=0, column=0)
max_long_label = tk.Label(aoi_frame, text="Maximum Longitude").grid(row=0, column=1)
min_lat_label = tk.Label(aoi_frame, text="Minimum Latitude").grid(row=0, column=2)
max_lat_label = tk.Label(aoi_frame, text="Maximum Latitude").grid(row=0, column=3)

min_long_entry = tk.Entry(aoi_frame)
max_long_entry = tk.Entry(aoi_frame)
min_lat_entry = tk.Entry(aoi_frame)
max_lat_entry = tk.Entry(aoi_frame)
min_long_entry.grid(row=1, column=0)
max_long_entry.grid(row=1, column=1)
min_lat_entry.grid(row=1, column=2)
max_lat_entry.grid(row=1, column=3)


for widget in aoi_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

# Satellite
sat_frame =tk.LabelFrame(main_tab, text="Satellite")
sat_frame.grid(row=2, column=0, sticky="news",padx=20, pady=10)

sat_combobox = ttk.Combobox(sat_frame, values=["Landsat 4", "Landsat 5", "Landsat 7", "Landsat 8", "Landsat 9", "Sentinel 1", "Sentinel 2"])
sat_combobox.grid(row=0, column=2)

for widget in sat_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

# Image processing
proc_frame = tk.LabelFrame(main_tab, text="Image processing")
proc_frame.grid(row=3, column=0, sticky="news", padx=20, pady=10)

reducer_label = tk.Label(proc_frame, text="Reducer").grid(row=0, column=0)
ccover_label = tk.Label(proc_frame, text="Maximum cloud cover %").grid(row=0, column=1)
cmask_label = tk.Label(proc_frame, text="Cloud mask").grid(row=0, column=2)

reducer_entry = ttk.Combobox(proc_frame, values=["Median", "Mosaic"])
reducer_entry.grid(row=1, column=0)
ccover_entry = tk.Spinbox(proc_frame, from_=0, to=100)
ccover_entry.grid(row=1, column=1)

cmask_flag = tk.IntVar(value=0)
cmask_check = tk.Checkbutton(proc_frame, variable=cmask_flag, onvalue=1, offvalue=0).grid(row=1, column=2)

for widget in proc_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

# Dates
dates_frame = tk.LabelFrame(main_tab, text="Dates")
dates_frame.grid(row=4, column=0, sticky="news", padx=20, pady=10)

sdate_label = tk.Label(dates_frame, text="Start Date").grid(row=0, column=0)
edate_label = tk.Label(dates_frame, text="End Date").grid(row=0, column=1)

sdate_entry = tk.Entry(dates_frame)
sdate_entry.grid(row=1, column=0)
edate_entry = tk.Entry(dates_frame)
edate_entry.grid(row=1, column=1)

deltp_label = tk.Label(dates_frame, text="Delta period").grid(row=2, column=0)
period_label = tk.Label(dates_frame, text="Period").grid(row=2, column=1)

deltp_entry = tk.Entry(dates_frame)
deltp_entry.grid(row=3, column=0)
period_entry = ttk.Combobox(dates_frame, values=["Day", "Week", "Month", "Year"])
period_entry.grid(row=3, column=1)

deltsp_label = tk.Label(dates_frame, text="Delta sub period").grid(row=4, column=0)
subperiod_label = tk.Label(dates_frame, text="Sub period").grid(row=4, column=1)
substart_label = tk.Label(dates_frame, text="Sub period start date").grid(row=4, column=2)

deltsp_entry = tk.Entry(dates_frame)
deltsp_entry.grid(row=5, column=0)
subperiod_entry = ttk.Combobox(dates_frame, values=["Day", "Week", "Month", "Year"])
subperiod_entry.grid(row=5, column=1)
substart_entry = tk.Entry(dates_frame)
substart_entry.grid(row=5, column=2)

for widget in dates_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

for widget in main_tab.winfo_children():
    widget.grid_configure(padx=10, pady=5)


#----------------------------------------
#             Masking tab
#----------------------------------------

# Masking
maskgen_frame = tk.LabelFrame(mask_tab, text="Generation")
maskgen_frame.grid(row=1, column=0, sticky="news", padx=20, pady=10)

band_label = tk.Label(maskgen_frame, text="Mask band").grid(row=0, column=0)
mupperb_label = tk.Label(maskgen_frame, text="Upper bound").grid(row=0, column=1)
mlowerb_label = tk.Label(maskgen_frame, text="Lower bound").grid(row=0, column=2)
rastf_label = tk.Label(maskgen_frame, text="Save as raster").grid(row=0, column=3)
vecf_label = tk.Label(maskgen_frame, text="Save as vector").grid(row=0, column=4)

mband_entry = tk.Entry(maskgen_frame)
mband_entry.grid(row=1, column=0)
mupperb_entry = tk.Entry(maskgen_frame)
mupperb_entry.grid(row=1, column=1)
mlowerb_entry = tk.Entry(maskgen_frame)
mlowerb_entry.grid(row=1, column=2)
rast_flag = tk.IntVar(value=0)
rastf_check = tk.Checkbutton(maskgen_frame, variable=rast_flag, onvalue=1, offvalue=0).grid(row=1, column=3)
vec_flag = tk.IntVar(value=0)
vecf_check = tk.Checkbutton(maskgen_frame, variable=vec_flag, onvalue=1, offvalue=0).grid(row=1, column=4)

for widget in maskgen_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

maskload_frame = tk.LabelFrame(mask_tab, text="Loading")
maskload_frame.grid(row=2, column=0, sticky="news", padx=20, pady=10)

mload_label = tk.Label(maskload_frame, text="Filepath").grid(row=0, column=0)
minv_label = tk.Label(maskload_frame, text="Invert").grid(row=0, column=1)
moverlap_label = tk.Label(maskload_frame, text="Overlap clip").grid(row=0, column=2)

mload_entry = tk.Entry(maskload_frame)
mload_entry.grid(row=1, column=0)
minv_flag = tk.IntVar(value=0)
minv_check = tk.Checkbutton(maskload_frame, variable=minv_flag, onvalue=1, offvalue=0).grid(row=1, column=1)
moverlap_entry = ttk.Combobox(maskload_frame, values=[1, 2])
moverlap_entry.grid(row=1, column=2)

for widget in maskload_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)
    
for widget in mask_tab.winfo_children():
    widget.grid_configure(padx=10, pady=5)


#----------------------------------------
#             Analysis tab
#----------------------------------------


# Indices
ind_frame = tk.LabelFrame(analysis_tab, text="Indices")
ind_frame.grid(row=0, column=0, sticky="news", padx=20, pady=10)

NDVI_label = tk.Label(ind_frame, text="NDVI").grid(row=0, column=0)
NDVI_flag = tk.IntVar(value=0)
NDVI_check = tk.Checkbutton(ind_frame, variable=NDVI_flag, onvalue=1, offvalue=0).grid(row=0, column=1,padx=10)

NDWI_label = tk.Label(ind_frame, text="NDWI").grid(row=1, column=0)
NDWI_flag = tk.IntVar(value=0)
NDWI_check = tk.Checkbutton(ind_frame, variable=NDWI_flag, onvalue=1, offvalue=0).grid(row=1, column=1)

NDMI_label = tk.Label(ind_frame, text="NDMI").grid(row=2, column=0)
NDMI_flag = tk.IntVar(value=0)
NDMI_check = tk.Checkbutton(ind_frame, variable=NDMI_flag, onvalue=1, offvalue=0).grid(row=2, column=1)

NDSI_label = tk.Label(ind_frame, text="NDSI").grid(row=3, column=0)
NDSI_flag = tk.IntVar(value=0)
NDSI_check = tk.Checkbutton(ind_frame, variable=NDSI_flag, onvalue=1, offvalue=0).grid(row=3, column=1)

NBR_label = tk.Label(ind_frame, text="NBR").grid(row=4, column=0)
NBR_flag = tk.IntVar(value=0)
NBR_check = tk.Checkbutton(ind_frame, variable=NBR_flag, onvalue=1, offvalue=0).grid(row=4, column=1)

NDRE_label = tk.Label(ind_frame, text="NDRE").grid(row=5, column=0)
NDRE_flag = tk.IntVar(value=0)
NDRE_check = tk.Checkbutton(ind_frame, variable=NDRE_flag, onvalue=1, offvalue=0).grid(row=5, column=1)

GNDVI_label = tk.Label(ind_frame, text="GNDVI").grid(row=6, column=0)
GNDVI_flag = tk.IntVar(value=0)
GNDVI_check = tk.Checkbutton(ind_frame, variable=GNDVI_flag, onvalue=1, offvalue=0).grid(row=6, column=1)

RECI_label = tk.Label(ind_frame, text="RECI").grid(row=7, column=0)
RECI_flag = tk.IntVar(value=0)
RECI_check = tk.Checkbutton(ind_frame, variable=RECI_flag, onvalue=1, offvalue=0).grid(row=7, column=1)

GCI_label = tk.Label(ind_frame, text="GCI").grid(row=8, column=0)
GCI_flag = tk.IntVar(value=0)
GCI_check = tk.Checkbutton(ind_frame, variable=GCI_flag, onvalue=1, offvalue=0).grid(row=8, column=1)

SIPI_label = tk.Label(ind_frame, text="SIPI").grid(row=9, column=0)
SIPI_flag = tk.IntVar(value=0)
SIPI_check = tk.Checkbutton(ind_frame, variable=SIPI_flag, onvalue=1, offvalue=0).grid(row=9, column=1)

EVI_label = tk.Label(ind_frame, text="EVI").grid(row=10, column=0)
EVI_flag = tk.IntVar(value=0)
EVI_check = tk.Checkbutton(ind_frame, variable=EVI_flag, onvalue=1, offvalue=0).grid(row=10, column=1)

SAVI_label = tk.Label(ind_frame, text="SAVI").grid(row=11, column=0)
SAVI_flag = tk.IntVar(value=0)
SAVI_check = tk.Checkbutton(ind_frame, variable=SAVI_flag, onvalue=1, offvalue=0).grid(row=11, column=1)

OSAVI_label = tk.Label(ind_frame, text="OSAVI").grid(row=12, column=0)
OSAVI_flag = tk.IntVar(value=0)
OSAVI_check = tk.Checkbutton(ind_frame, variable=OSAVI_flag, onvalue=1, offvalue=0).grid(row=12, column=1)

MSAVI_label = tk.Label(ind_frame, text="MSAVI").grid(row=13, column=0)
MSAVI_flag = tk.IntVar(value=0)
MSAVI_check = tk.Checkbutton(ind_frame, variable=MSAVI_flag, onvalue=1, offvalue=0).grid(row=13, column=1)

ARVI_label = tk.Label(ind_frame, text="ARVI").grid(row=14, column=0)
ARVI_flag = tk.IntVar(value=0)
ARVI_check = tk.Checkbutton(ind_frame, variable=ARVI_flag, onvalue=1, offvalue=0).grid(row=14, column=1)

VARI_label = tk.Label(ind_frame, text="VARI").grid(row=15, column=0)
VARI_flag = tk.IntVar(value=0)
VARI_check = tk.Checkbutton(ind_frame, variable=VARI_flag, onvalue=1, offvalue=0).grid(row=15, column=1)

ARI_label = tk.Label(ind_frame, text="ARI").grid(row=16, column=0)
ARI_flag = tk.IntVar(value=0)
ARI_check = tk.Checkbutton(ind_frame, variable=ARI_flag, onvalue=1, offvalue=0).grid(row=16, column=1)

MARI_label = tk.Label(ind_frame, text="MARI").grid(row=17, column=0)
MARI_flag = tk.IntVar(value=0)
MARI_check = tk.Checkbutton(ind_frame, variable=MARI_flag, onvalue=1, offvalue=0).grid(row=17, column=1)

BSI_label = tk.Label(ind_frame, text="BSI").grid(row=18, column=0)
BSI_flag = tk.IntVar(value=0)
BSI_check = tk.Checkbutton(ind_frame, variable=BSI_flag, onvalue=1, offvalue=0).grid(row=18, column=1)

PSRI_label = tk.Label(ind_frame, text="PSRI").grid(row=19, column=0)
PSRI_flag = tk.IntVar(value=0)
PSRI_check = tk.Checkbutton(ind_frame, variable=PSRI_flag, onvalue=1, offvalue=0).grid(row=19, column=1)

OSI_label = tk.Label(ind_frame, text="OSI").grid(row=20, column=0)
OSI_flag = tk.IntVar(value=0)
OSI_check = tk.Checkbutton(ind_frame, variable=OSI_flag, onvalue=1, offvalue=0).grid(row=20, column=1)

# RGB Composites
comp_frame = tk.LabelFrame(analysis_tab, text="RGB Composites")
comp_frame.grid(row=0, column=1, sticky="news", padx=20, pady=10)

tcol_label = tk.Label(comp_frame, text="True colour").grid(row=0, column=0)
RGB_flag = tk.IntVar(value=0)
tcol_check = tk.Checkbutton(comp_frame, variable=RGB_flag, onvalue=1, offvalue=0).grid(row=0, column=1,padx=10)

rveg1_label = tk.Label(comp_frame, text="Red vegetation 1").grid(row=1, column=0)
RedVegetation1_flag = tk.IntVar(value=0)
rveg1_check = tk.Checkbutton(comp_frame, variable=RedVegetation1_flag, onvalue=1, offvalue=0).grid(row=1, column=1)

rveg2_label = tk.Label(comp_frame, text="Red vegetation 2").grid(row=2, column=0)
RedVegetation2_flag = tk.IntVar(value=0)
rveg2_check = tk.Checkbutton(comp_frame, variable=RedVegetation2_flag, onvalue=1, offvalue=0).grid(row=2, column=1)

SWIR1_label = tk.Label(comp_frame, text="SWIR composite 1").grid(row=3, column=0)
SWIRComposite1_flag = tk.IntVar(value=0)
SWIR1_check = tk.Checkbutton(comp_frame, variable=SWIRComposite1_flag, onvalue=1, offvalue=0).grid(row=3, column=1)

SWIR2_label = tk.Label(comp_frame, text="SWIR composite 2").grid(row=4, column=0)
SWIRComposite2_flag = tk.IntVar(value=0)
SWIR2_check = tk.Checkbutton(comp_frame, variable=SWIRComposite2_flag, onvalue=1, offvalue=0).grid(row=4, column=1)

Agri_label = tk.Label(comp_frame, text="Agriculture").grid(row=5, column=0)
Agriculture_flag = tk.IntVar(value=0)
Agri_check = tk.Checkbutton(comp_frame, variable=Agriculture_flag, onvalue=1, offvalue=0).grid(row=5, column=1)

bathy1_label = tk.Label(comp_frame, text="Bathymetric 1").grid(row=6, column=0)
Bathymetric1_flag = tk.IntVar(value=0)
bathy1_check = tk.Checkbutton(comp_frame, variable=Bathymetric1_flag, onvalue=1, offvalue=0).grid(row=6, column=1)

bathy2_label = tk.Label(comp_frame, text="Bathymetric 2").grid(row=7, column=0)
Bathymetric2_flag = tk.IntVar(value=0)
bathy2_check = tk.Checkbutton(comp_frame, variable=Bathymetric2_flag, onvalue=1, offvalue=0).grid(row=7, column=1)

bathy3_label = tk.Label(comp_frame, text="Bathymetric 3").grid(row=8, column=0)
Bathymetric3_flag = tk.IntVar(value=0)
bathy3_check = tk.Checkbutton(comp_frame, variable=Bathymetric3_flag, onvalue=1, offvalue=0).grid(row=8, column=1)

bathy4_label = tk.Label(comp_frame, text="Bathymetric 4").grid(row=9, column=0)
Bathymetric4_flag = tk.IntVar(value=0)
bathy4_check = tk.Checkbutton(comp_frame, variable=Bathymetric4_flag, onvalue=1, offvalue=0).grid(row=9, column=1)

bathy5_label = tk.Label(comp_frame, text="Bathymetric 5").grid(row=10, column=0)
Bathymetric5_flag = tk.IntVar(value=0)
bathy5_check = tk.Checkbutton(comp_frame, variable=Bathymetric5_flag, onvalue=1, offvalue=0).grid(row=10, column=1)

bathy6_label = tk.Label(comp_frame, text="Bathymetric 6").grid(row=11, column=0)
Bathymetric6_flag = tk.IntVar(value=0)
bathy6_check = tk.Checkbutton(comp_frame, variable=Bathymetric6_flag, onvalue=1, offvalue=0).grid(row=11, column=1)

urban_label = tk.Label(comp_frame, text="Urban").grid(row=12, column=0)
Urban_flag = tk.IntVar(value=0)
urban_check = tk.Checkbutton(comp_frame, variable=Urban_flag, onvalue=1, offvalue=0).grid(row=12, column=1)

Neon_label = tk.Label(comp_frame, text="Neon").grid(row=13, column=0)
Neon_flag = tk.IntVar(value=0)
Neon_check = tk.Checkbutton(comp_frame, variable=Neon_flag, onvalue=1, offvalue=0).grid(row=13, column=1)

# Quantitative analyses
quant_frame = tk.LabelFrame(analysis_tab, text="Quantitative analyses")
quant_frame.grid(row=0, column=2, sticky="news", padx=20, pady=10)

areas_label = tk.Label(quant_frame, text="Areas").grid(row=0, column=0)
areas_entry = tk.Entry(quant_frame)
areas_entry.grid(row=1, column=0)

band_label = tk.Label(quant_frame, text="Band name").grid(row=2, column=0)
band_entry = tk.Entry(quant_frame)
band_entry.grid(row=3, column=0)

mean_label = tk.Label(quant_frame, text="Area mean").grid(row=4, column=0)
AreaMean_flag = tk.IntVar(value=0)
mean_check = tk.Checkbutton(quant_frame, variable=AreaMean_flag, onvalue=1, offvalue=0).grid(row=4, column=1)

stddev_label = tk.Label(quant_frame, text="Area standard deviation").grid(row=5, column=0)
AreaStdDev_flag = tk.IntVar(value=0)
stddev_check = tk.Checkbutton(quant_frame, variable=AreaStdDev_flag, onvalue=1, offvalue=0).grid(row=5, column=1)

thresh_label = tk.Label(quant_frame, text="Area threshold").grid(row=6, column=0)
AreaThreshold_flag = tk.IntVar(value=0)
thresh_check = tk.Checkbutton(quant_frame, variable=AreaThreshold_flag, onvalue=1, offvalue=0).grid(row=6, column=1)

pval_label = tk.Label(quant_frame, text="Point value").grid(row=7, column=0)
PointValue_flag = tk.IntVar(value=0)
pval_check = tk.Checkbutton(quant_frame, variable=PointValue_flag, onvalue=1, offvalue=0).grid(row=7, column=1)

qaupperb_label = tk.Label(quant_frame, text="Upper bound").grid(row=8, column=0)
qalowerb_label = tk.Label(quant_frame, text="Lower bound").grid(row=10, column=0)
qaupperb_entry = tk.Entry(quant_frame)
qaupperb_entry.grid(row=9, column=0)
qalowerb_entry = tk.Entry(quant_frame)
qalowerb_entry.grid(row=11, column=0)

qasave_label = tk.Label(quant_frame, text="Show areas").grid(row=12, column=0)
qasave_flag = tk.IntVar(value=0)
qasave_check = tk.Checkbutton(quant_frame, variable=qasave_flag, onvalue=1, offvalue=0).grid(row=12, column=1)

for widget in analysis_tab.winfo_children():
    widget.grid_configure(padx=10, pady=5)

#Reading in all the data
def enter_data():
    #Getting all the input data
    run_name, imsave_name, masksave_name, plotsave_name = run_name_entry.get(), imsave_name_entry.get(), masksave_name_entry.get(), plotsave_name_entry.get()
    min_long, max_long, min_lat, max_lat = min_long_entry.get(), max_long_entry.get(), min_lat_entry.get(), max_lat_entry.get()
    satellite = sat_combobox.get()
    reducer, ccover, cmask = reducer_entry.get(), ccover_entry.get(), cmask_flag.get()
    sdate, edate, deltp, period, deltsp, subperiod, substart = sdate_entry.get(), edate_entry.get(), deltp_entry.get(), period_entry.get(), deltsp_entry.get(), subperiod_entry.get(), substart_entry.get()
    
    mband, mupper, mlower, rast_f, vec_f = mband_entry.get(), mupperb_entry.get(), mlowerb_entry.get(), rast_flag.get(), vec_flag.get() 
    mload, minv_f, moverlap = mload_entry.get(), minv_flag.get(), moverlap_entry.get()

    #The analysis lists require slightly different handling to create lists from a series of flags
    ans_inds_bin = [globals()[an+'_flag'].get() for an in Analysis_lib]
    ans_inds = [i for i, x in enumerate(ans_inds_bin) if x == 1]
    ans_list = [Analysis_lib[ind] for ind in ans_inds]

    areas, qband, qaupper, qalower, qasave_f = areas_entry.get(), band_entry.get(), qaupperb_entry.get(), qalowerb_entry.get(), qasave_flag.get() 
    qans_inds_bin = [globals()[an+'_flag'].get() for an in Quant_lib]
    qans_inds = [i for i, x in enumerate(qans_inds_bin) if x == 1]
    qans_list = [Quant_lib[ind] for ind in qans_inds]
  
    #The run card giant string itself
    yaml_str = """\
#---------------------------------------------------------------------
#                  Google Earth Engine Package                       
#                                                                    
#  This file is used to set the parameters of the run.               
#                                                                    
#  Some notation/conventions:                                        
#                                                                    
#   Lines starting with a '# ' are info or comments                  
#                                                                    
#   Other lines follow:   variable: value     # comment  
#   All flags follow 0/1 for off/on           
#                                                   
#---------------------------------------------------------------------
        
                                                        
#---------------------------------------------------------------------
# Information for saving files                                       
#---------------------------------------------------------------------
Names:
  Name: """+run_name+"""                   # Main name of the run, prefixes all filenames 
  ImSaveDir: """+imsave_name+"""              # Google Drive folder images will be saved to
  MaskSaveDir: """+masksave_name+"""            # Google Drive folder any masks will be saved to
  PlotDir: """+plotsave_name+"""                # Local folder where plots will be saved
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Area of interest                               
#---------------------------------------------------------------------
min_long: """+min_long+"""         # Minimum longitude of the AOI 
max_long: """+max_long+"""         # Maximum longitude of the AOI 
min_lat: """+min_lat+"""          # Minimum latitude of the AOI 
max_lat: """+max_lat+"""          # Maximum latitude of the AOI 
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Satellite                               
#---------------------------------------------------------------------
Satellite: """+satellite+"""        # Name of the satellite to use data from
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Image processing                               
#---------------------------------------------------------------------
Processing:
  reducer: """+reducer+"""              # Operation by which to reduce individual collections to single images
  cloud_cover: """+str(ccover)+"""          # Cloud cover filter level
  cloud_mask: """+str(cmask)+"""           # Cloud mask flag
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Date information - format dates as 'YYYY-MM-dd'                              
#---------------------------------------------------------------------
Dates:
  start_date: """+str(sdate)+"""              # Start of the overall period
  end_date: """+str(edate)+"""                # End of the overall period
  delta_period: """+str(deltp)+"""            # Units of period time to sample
  period: """+period+"""                  # Time unit for the period
  delta_sub: """+str(deltsp)+"""               # Units of sub-period time to sample
  subperiod: """+subperiod+"""               # Time unit for the sub-period
  sub_start: """+substart+"""               # Start of the sub-period to sample 
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Masking criteria                        
#---------------------------------------------------------------------
Masking:
  mask_bandName: """+mband+"""        # Name of the band used to calculate the mask
  upper_bound: """+str(mupper)+"""          # Upper bound of the mask
  lower_bound: """+str(mlower)+"""          # Lower bound of the mask
  save_raster: """+str(rast_f)+"""          # Flag for saving the mask as a raster image
  save_vector: """+str(vec_f)+"""          # Flag for saving the mask as a vector image (.shp)
  load_name: """+mload+"""            # Filepath of an Earth Engine shapefile asset to load
  load_invert: """+str(minv_f)+"""          # Flag for inversion of a loaded mask
  load_overlap: """+str(moverlap)+"""         # [1,2] - Length of overlap regions to clip from a loaded mask     
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Analytic information - which indices to calculate and analyses to perform                          
#---------------------------------------------------------------------
Analyses:
  indices: """+str(ans_list)+"""               # List of indices to calculate over the whole image
  areas: """+str(areas)+"""                 # List of coordinates lists, either points or geometries, to be analysed
  bandName: """+qband+"""              # Band on which to perform the quantitative analyses
  quant_analyses: """+str(qans_list)+"""        # Quantitative analyses to be performed
  upper_bound: """+str(qaupper)+"""           # Upper bound for a threshold analysis
  lower_bound: """+str(qalower)+"""           # Lower bound for a threshold analysis
  save_flag: """+str(qasave_f)+"""             # Flag for saving an RGB image with analysed areas highlighted
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Visualisations                    
#---------------------------------------------------------------------
Visualisations:

#---------------------------------------------------------------------
"""
    print(sdate)
    print(str(sdate))

    #Writing to file
    with open('RunCards/'+run_name+'.yml', "w",) as yaml_file:
        yaml_file.write(yaml_str)

    #Checking the user inputs
    print('Checking user inputs...')
    names_l, aoi_l, sat_l, process_dict_l, dates_dict_l, mask_dict_l, analyses_dict_l, vis_dict_l = RunCardReader(run_name+'.yml')

    #If all valid, close the GUI and pass on the information
    if names_l is not None:
        root.destroy()
        global names 
        global aoi
        global sat
        global process_dict
        global dates_dict
        global mask_dict
        global analyses_dict
        global vis_dict
        names, aoi, sat, process_dict, dates_dict, mask_dict, analyses_dict, vis_dict = names_l, aoi_l, sat_l, process_dict_l, dates_dict_l, mask_dict_l, analyses_dict_l, vis_dict_l

        return 
    
    

# Data entry button
button1 = ttk.Button(main_tab, text="Enter data", command= enter_data).grid(row=60, column=0, sticky="news", padx=20, pady=10)
button2 = ttk.Button(mask_tab, text="Enter data", command= enter_data).grid(row=60, column=0, sticky="news", padx=20, pady=10)
button3 = ttk.Button(analysis_tab, text="Enter data", command= enter_data).grid(row=60, column=0, sticky="news", padx=20, pady=10)

#Generate the GUI
root.mainloop()  