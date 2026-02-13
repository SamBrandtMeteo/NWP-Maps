##########################
# SECTION 1: USER INPUTS #
##########################

# Time (in UTC)
# The selected time must be contained within the latest GFS run.
year = 2026
month = 2
day = 13
hour = 12

# Edges of domain (in deg north of latitude and deg east of longitude)
north = 56
south = 21
east = -65
west = -130

# Cloud fraction threshold
CF = 25

########################
# SECTION 2: LIBRARIES #
########################

# The following Python libraries are required:
# numpy
# matplotlib
# siphon
# xarray
# scipy
# cartopy

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from siphon.catalog import TDSCatalog
from xarray import open_dataset
from xarray.backends import NetCDF4DataStore
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#############################
# SECTION 5: DATA RETRIEVAL #
#############################

# Define location of latest GFS run on THREDDS
THREDDS = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/'
url = THREDDS+'GFS/Global_0p25deg/latest.xml'

# Define location of the data
latest_model = TDSCatalog(url)
latest_ds = list(latest_model.datasets.values())[0]
ncss = latest_ds.subset()

# Create a datetime object to specify the output time that you want
valid = datetime(year,month,day,hour)

# Establish a query for the data
query = ncss.query()

# Trim data to location/time of interest
query.lonlat_box(north=north,south=south,east=east,west=west).time(valid)

# Specify that output needs to be in netcdf format
query.accept('netcdf4')

# Specify the variables that you want
var1 = 'Total_cloud_cover_isobaric'
var2 = 'Geopotential_height_isobaric'
var3 = 'Geopotential_height_surface'
query.variables(var1,var2,var3)

# Retrieve the data using the info from the query
data = ncss.get_data(query)

# Open the dataset
data = open_dataset(NetCDF4DataStore(data))

# Specify coordinate grids
lat = np.flip(np.arange(south,north+0.25,0.25))
lon = np.arange(west,east+0.25,0.25)
lon,lat = np.meshgrid(lon,lat)

# Import relevant variables
surface_height = np.array(data['Geopotential_height_surface'])[0,:,:]
cloud_fraction = np.flip(np.array(data['Total_cloud_cover_isobaric'])[0,:,:,:],axis=0)
isobaric_height = np.flip(np.array(data['Geopotential_height_isobaric'])[0,-22:,:,:],axis=0)

###########################
# SECTION 7: CALCULATIONS #
###########################

# Subtract surface height
height = isobaric_height-surface_height

# Calculate indices where cloud fraction exceeds threshold
ceiling_indices = np.argmax(cloud_fraction>=CF,axis=0)

# Apply indices to the height grid
ceiling_height = np.take_along_axis(height,ceiling_indices[None, :, :],axis=0)[0]

# Mask grid points with no cloud cover
ceiling_height = np.ma.masked_where(np.logical_and(ceiling_indices==0,cloud_fraction[0,:,:]<CF),ceiling_height)
#ceiling_height=np.ma.masked_where(ceiling_indices==0,ceiling_height)

#######################
# SECTION 7: PLOTTING #
#######################

# Create geographic axis
proj_kwargs={'projection': ccrs.PlateCarree(),'adjustable': 'box'}
fig,ax=plt.subplots(subplot_kw=proj_kwargs,dpi=500)

# Add country borders
ax.add_feature(cfeature.COASTLINE.with_scale('50m'),edgecolor='gray',lw=0.5)

# Add state/territory borders
ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray',lw=0.5)

# Filled contours of 500 mb geostrophic vorticity
ct_vort=ax.pcolormesh(lon,lat,ceiling_height/1000,cmap='plasma',vmin=0,vmax=12)

# Add colorbar
cbar=plt.colorbar(ct_vort,ax=ax,orientation='horizontal',fraction=0.04,pad=0.02)

# Set colorbar tick label size
cbar.ax.tick_params(labelsize=6)

# Add colorbar label
cbar.ax.set_xlabel('Ceiling Height ($km$)',fontsize=6)

# Add figure title
tit1='Ceiling Height Where Cloud Fraction Exceeds '+str(CF)+'%\n'
init1='GFS Initialized '+ncss.metadata.time_span['begin'][0:10]
init2=' '+ncss.metadata.time_span['begin'][11:13]+'z, '
valid='Valid '+str(valid)[0:13]+'z'
ax.set_title(tit1+init1+init2+valid,fontsize=8)
sig='https://github.com/SamBrandtMeteo'
ax.text(0.99,0.01,sig,transform=ax.transAxes,ha='right',va='bottom',fontsize=4,color='black')

plt.show()
