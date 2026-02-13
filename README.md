
NWP-Maps is a collection of Python scripts that plot interesting things with NWP data

-------------------------------------------------------------------------------------

IDGZ.py

IDGZ is a parameter I came up with called the "Integrated Dendritic Growth Zone". It
is equal to relative humidity multiplied by omega (vertical velocity in pressure
coordinates), integrated from the bottom of the DGZ at -12 C to the top of the DGZ
at -17 C. It is meant to be used as a quick way to determine where favorabe 
conditions for dendrite growth exist. IDGZ is maximized where the depth of the DGZ,
ascent within the DGZ, and/or relative humidity within the DGZ are largest.

-------------------------------------------------------------------------------------

BaroclinicWave.py

Baroclinic instability allows synoptic waves to couple with low pressure systems and 
amplify in the presence of a horizontal temperature gradient. The process involves a
positive feedback loop between low-level pressure tendiencies from differential 
cyclonic vortictiy advection (DCVA) and upper level height tendencies from low-level 
temperature advection. BaroclinicWave.py plots the 500 mb heights and geostrophic 
vorticity under the assumption that 500 mb is near the level of non-divergence (LND),
and thus vertical motion from DCVA is approximately maximized. It would technically
be more accurate to explicitly show the vertical gradient of vorticity advection, but
the 500 mb CVA approximation tends to work quite well. The second part of these maps
is 850 mb geostrophic temperature advection. Since temperature advection tends to be
maximized near the surface, the quasi-geostrophic height tendency says that mid-level
heights should fall above low-level cold advection and vice versa. To identify the
release of baroclinic instability, look for couplets of warm and cold advection in
the vicinity of CVA, which implies wave amplification.

-------------------------------------------------------------------------------------

CloudCeiling.py

The cloud ceiling, or the lowest height above the ground where cloud cover first
occurs, is a useful parameter for sunrise/sunset forecasting. Larger ceilings allow
more sunlight to illuminate the cloud base from below, producing more vibrant
displays.
