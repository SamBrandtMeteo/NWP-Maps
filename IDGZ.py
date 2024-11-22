##########################
# SECTION 1: USER INPUTS #
##########################

# Model Select (options are 'RAP' and 'GFS')
# RAP: 13 km grid spacing, good for real-time mesoanalysis
# GFS: 0.25 degree grid spacing, good for medium range forecasting
model='RAP'

# Time (in UTC)
# The selected time must be contained within the latest model run.
year=2024
month=11
day=22
hour=18

# Edges of Domain (in degrees north of latitude and degrees east of longitude)
# GFS covers entire globe, RAP only covers CONUS
north=45
south=35
east=-70
west=-85

# Option to save the map to a specified location on your computer.
# If set to True, define the location parameter below.
savefig=False

if savefig==True:

  # Specify where on your computer to save the file
  location='.../your_directory_here/'

  # Filename
  date=str(year)+str(month).zfill(2)+str(day).zfill(2)
  filename='IGDZ_'+model+'_'+date+'_'+str(hour).zfill(2)+'z.png'

########################
# SECTION 2: LIBRARIES #
########################

# The following Python libraries are required:
# numpy
# matplotlib
# siphon
# xarray
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
# SECTION 3: DATA RETRIEVAL #
#############################

# Define location of latest model run on THREDDS
THREDDS='https://thredds.ucar.edu/thredds/catalog/grib/NCEP/'
if model=='GFS':
  url=THREDDS+'GFS/Global_0p25deg/latest.xml'
elif model=='RAP':
  url=THREDDS+'RAP/CONUS_13km/latest.xml'

# Define location of the data
latest_model=TDSCatalog(url)
latest_ds=list(latest_model.datasets.values())[0]
ncss=latest_ds.subset()

# Create a datetime object to specify the output time that you want
valid=datetime(year,month,day,hour)

# Establish a query for the data
query=ncss.query()

# Trim data to location/time of interest
query.lonlat_box(north=north,south=south,east=east,west=west).time(valid)

# Include irregular lat/lon grid for RAP case
if model=='RAP':
  query.add_lonlat(value=True)

# Specify that output needs to be in netcdf format
query.accept('netcdf4')

# Specify the variables that you want
var1='Categorical_Snow_surface'
var2='Temperature_isobaric'
var3='Vertical_velocity_pressure_isobaric'
var4='Relative_humidity_isobaric'
query.variables(var1,var2,var3,var4)

# Retrieve the data using the info from the query
data=ncss.get_data(query)

# Open the dataset
data=open_dataset(NetCDF4DataStore(data))

# Specify coordinate grids based on model choice
if model=='GFS':
  metadata_isobaric=ncss.metadata.axes['isobaric']['attributes'][2]['values']
  plevs=np.array(list(map(float,metadata_isobaric)))
  lat=np.flip(np.arange(south,north+0.25,0.25))
  lon=np.arange(west,east+0.25,0.25)
  lon,lat=np.meshgrid(lon,lat)
else:
  plevs=np.arange(100,1025,25)*100
  lon=np.array(data['lon'])
  lat=np.array(data['lat'])

# Import snow ID grid
snow=np.array(data['Categorical_Snow_surface'])[0,:,:]

# Import temperature grid
T=np.array(data['Temperature_isobaric'])[0,plevs>=250,:,:]

# Import vertical velocity grid
Omega=np.array(data['Vertical_velocity_pressure_isobaric'])[0,plevs>=250,:,:]

# Import relative humidity grid
RH=np.array(data['Relative_humidity_isobaric'])[0,plevs>=250,:,:]/100

# Trim pressure levels to prevent double-counting stratospheric contributions 
# to IDGZ. This is really only a possible issue with GFS isobaric data, but
# the trimming is also done with the RAP data to streamline the code.
plevs=plevs[plevs>=250]

###########################
# SECTION 4: CALCULATIONS #
###########################

# Create empty IDGZ grid
IDGZ=np.zeros(np.shape(T[0,:,:]))

# Loop through every lat/lon
for i in range(0,len(T[0,:,0])):
    for j in range(0,len(T[0,0,:])):

        # Find the part of the column that is entirely within the DGZ
        indices=np.where(np.logical_and(T[:,i,j]>256.15,T[:,i,j]<261.15))[0]

        # Add the DGZ boundary conditions
        T_DGZ=np.array([256.15])
        T_DGZ=np.append(T_DGZ,T[indices,i,j])
        T_DGZ=np.append(T_DGZ,261.15)

        # Interpolate other variables to DGZ boundaries
        # Results are flipped so the integration goes from bottom to top
        p_DGZ=np.flip(np.interp(T_DGZ,T[:,i,j],plevs))
        Omega_I=np.flip(np.interp(p_DGZ,plevs,Omega[:,i,j]))
        RH_I=np.flip(np.interp(p_DGZ,plevs,RH[:,i,j]))

        # Calculate IDGZ using trapezoidal integration
        IDGZ[i,j]=np.trapz(RH_I*Omega_I,p_DGZ)

# Mask IDGZ where the model says it is not snowing
IDGZ=np.ma.masked_where(snow==0,IDGZ)

#######################
# SECTION 5: PLOTTING #
#######################

# Create geographic axis
proj_kwargs={'projection': ccrs.PlateCarree(),'adjustable': 'box'}
fig,ax=plt.subplots(subplot_kw=proj_kwargs,dpi=500)

# Add country borders
ax.add_feature(cfeature.COASTLINE.with_scale('50m'),edgecolor='gray',lw=0.5)

# Add state/territory borders
ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray',lw=0.5)

# Filled contour plot of IDGZ
ctf=ax.contourf(lon,lat,IDGZ/1000,np.arange(0,33,3),cmap='BuPu',extend='max')

# Add colorbar
cbar = plt.colorbar(ctf,ax=ax,orientation='horizontal',fraction=0.04,pad=0.02)

# Set colorbar tick label size
cbar.ax.tick_params(labelsize=6)

# Add colorbar label
cbar.ax.set_xlabel('IDGZ'+' ($10^{-3}$ $Pa^{2}$ $s^{-1}$)',fontsize=6)

# Add figure title
var='Integrated Dendritic Growth Zone\n'
init1=' Initialized '+ncss.metadata.time_span['begin'][0:10]
init2=' '+ncss.metadata.time_span['begin'][11:13]+'z, '
valid='Valid '+str(valid)[0:13]+'z'
ax.set_title(var+model+init1+init2+valid,fontsize=8)
sig='https://github.com/SamBrandtMeteo'
ax.text(0.99,0.01,sig,transform=ax.transAxes,ha='right',va='bottom',fontsize=4)

# If selected, save map
if savefig==True:
  plt.savefig(location+filename,bbox_inches='tight')
