##########################
# SECTION 1: USER INPUTS #
##########################

# Time (in UTC)
# The selected time must be contained within the latest GFS run.
year=2024
month=11
day=24
hour=0

# Edges of Domain (in degrees north of latitude and degrees east of longitude)
north=55
south=25
east=-70
west=-100

# Option to save the map to a specified location on your computer.
# If set to True, define the location parameter below.
savefig=False

if savefig==True:

  # Specify where on your computer to save the file
  location='.../your_directory_here/'

  # Filename
  date=str(year)+str(month).zfill(2)+str(day).zfill(2)
  filename='BaroclinicWave_'+model+'_'+date+'_'+str(hour).zfill(2)+'z.png'

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
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#################################
# SECTION 3: PHYSICAL CONSTANTS #
#################################

g=9.81 # Acceleration due to gravity in m/s^2
Omega=7.2921*10**-5 # Rotation rate of Earth in rad/s

########################
# SECTION 4: FUNCTIONS #
########################

def partial(lat,lon,field,wrt):
    # Calculates horizontal gradients on a lat/lon grid
    gradient=np.zeros(np.shape(field))
    if wrt=='x':
        upper=field[:,2::]
        lower=field[:,0:-2]
        dx=111200*np.cos(lat[:,2::]*(np.pi/180))*(lon[0,1]-lon[0,0])
        grad=(upper-lower)/(2*dx)
        gradient[:,1:-1]=grad
        gradient[:,0]=grad[:,0]
        gradient[:,-1]=grad[:,-1]
    if wrt=='y':
        upper=field[2::,:]
        lower=field[0:-2,:]
        dy=111200*(lat[1,0]-lat[0,0])
        grad=(upper-lower)/(2*dy)
        gradient[1:-1,:]=grad
        gradient[0,:]=grad[0,:]
        gradient[-1,:]=grad[-1,:]
    return gradient

def laplacian(lat,lon,field):
    # Calculates the horizontal laplacian on a lat/lon grid
    d2xdx2=partial(lat,lon,partial(lat,lon,field,'x'),'x')
    d2ydy2=partial(lat,lon,partial(lat,lon,field,'y'),'y')
    return d2xdx2+d2ydy2

def advection(lat,lon,u,v,field):
    # Calculates horizontal advection of a field on a lat/lon grid
    return u*partial(lat,lon,field,'x')+v*partial(lat,lon,field,'y')

#############################
# SECTION 5: DATA RETRIEVAL #
#############################

# Define location of latest GFS run on THREDDS
THREDDS='https://thredds.ucar.edu/thredds/catalog/grib/NCEP/'
url=THREDDS+'GFS/Global_onedeg/latest.xml'

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

# Specify that output needs to be in netcdf format
query.accept('netcdf4')

# Specify the variables that you want
var1='Geopotential_height_isobaric'
var2='Temperature_isobaric'
query.variables(var1,var2)

# Retrieve the data using the info from the query
data=ncss.get_data(query)

# Open the dataset
data=open_dataset(NetCDF4DataStore(data))

# Specify coordinate grids
metadata_isobaric=ncss.metadata.axes['isobaric']['attributes'][2]['values']
plevs=np.array(list(map(float,metadata_isobaric)))
lat=np.flip(np.arange(south,north+1))
lon=np.arange(west,east+1)
lon,lat=np.meshgrid(lon,lat)

# Define indices for the pressure levels of interest
index_850=np.where(plevs==850*100)[0][0]
index_500=np.where(plevs==500*100)[0][0]

# Import 500 mb geopotential height grid
z_500=np.array(data['Geopotential_height_isobaric'][0,index_500,:,:])

# Import 850 mb geopotential height grid
z_850=np.array(data['Geopotential_height_isobaric'][0,index_850,:,:])

# Import 850 mb temperature grid
T_850=np.array(data['Temperature_isobaric'][0,index_850,:,:])

###########################
# SECTION 6: CALCULATIONS #
###########################

# Calculate planetary vertical vorticity 
f=2*Omega*np.sin(lat*(np.pi/180))

# Calculate 850 mb geostrophic wind
ug_850=-(g/f)*partial(lat,lon,z_850,'y')
vg_850=(g/f)*partial(lat,lon,z_850,'x')

# Calculate 500 mb geostrophic vertical vorticity
vort_500=(g/f)*laplacian(lat,lon,z_500)

# Apply smoothing and convert to dam
gaussian_filter(z_500/10,sigma=1)

# Calculate 850 geostrophic temperature advection
adv_850=advection(lat,lon,ug_850,vg_850,T_850)

# Apply smoothing and convert to K/hr
adv_850=-gaussian_filter(adv_850,sigma=1)*3600

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
ct_vort=ax.contourf(lon,lat,gaussian_filter(vort_500,sigma=1)*10**5,
                    np.arange(0,65,5),cmap='plasma_r',extend='max')

# Contours of 500 mb heights
ct_500=ax.contour(lon,lat,gaussian_filter(z_500/10,sigma=1),
                  np.arange(0,10000,6),linewidths=0.5,colors='black',zorder=3)

# Contours of 850 mb geostrophic cold advection
ct_adv_cold=ax.contour(lon,lat,adv_850,np.arange(-10,0,0.5),
                       colors='deepskyblue',linewidths=0.75,zorder=4,
                       linestyles='solid')

# Contours of 850 mb geostrophic warm advection
ct_adv_warm=ax.contour(lon,lat,adv_850,np.arange(0.5,10.5,0.5),colors='red',
                       linewidths=0.75,zorder=4)

# Add contour labels
ax.clabel(ct_500,ct_500.levels,fontsize=4)
ax.clabel(ct_adv_cold,ct_adv_cold.levels,fontsize=4)
ax.clabel(ct_adv_warm,ct_adv_warm.levels,fontsize=4)

# Add colorbar
cbar=plt.colorbar(ct_vort,ax=ax,orientation='horizontal',fraction=0.04,pad=0.02)

# Set colorbar tick label size
cbar.ax.tick_params(labelsize=6)

# Add colorbar label
cbar.ax.set_xlabel('Vorticity ($10^{-5}$ $s^{-1}$)',fontsize=6)

# Add figure title
tit1='500 mb Heights (Black, $dam$) '
tit2='& Geostrophic Vorticity (Fill, $10^{-5}$ $s^{-1}$)\n'
tit3='850 mb Geostrophic Temperature Advection (Red and Blue, $K$ $hr^{-1}$)\n'
init1=' Initialized '+ncss.metadata.time_span['begin'][0:10]
init2=' '+ncss.metadata.time_span['begin'][11:13]+'z, '
valid='Valid '+str(valid)[0:13]+'z'
ax.set_title(tit1+tit2+tit3+init1+init2+valid,fontsize=8)
sig='https://github.com/SamBrandtMeteo'
ax.text(0.99,0.01,sig,transform=ax.transAxes,ha='right',va='bottom',fontsize=4)

# If selected, save map
if savefig==True:
  plt.savefig(location+filename,bbox_inches='tight')
