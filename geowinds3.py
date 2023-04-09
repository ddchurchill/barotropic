import cmath # used for complex numbers - storing wind vectors
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.interpolate
from velocity import Velocity as vel
from trajectory import Trajectory 
from traject import euler
from traject import huen

def get_zeta(z, x, y):
    """
    get_zeta - compute relative vorticity from geopotential height, z, with input longitudes, x,

    and input latitudes, y, 
    Use manual centered differences.
    Assuming grid spacing of 1 degree
    1st dimension is 60 of latitude, 2nd dimension is longitude
    This produces vorticity that agrees with the code using gradient calls to within 2%.
    
    """
    d2zdx2 = np.zeros((60, 140))
    d2zdy2 = np.zeros((60,140))
    for i in range(2,138):
        for j in range(2, 59) :
            dlon = (x[j,i+1] - x[j,i-1])/2.
#            print("dlon: ", dlon)
            dx = earth_radius * np.cos(np.radians(y[j,i])) * np.radians( dlon)
            d2zdx2[j,i] = (z[j,i+1] - 2*z[j,i] + z[j,i-1])/dx**2
 
            dlat = (y[j+1,i] - y[j-1,i])/2.
#            print("dlat:", dlat)
            dy = earth_radius * np.radians(dlat)
            d2zdy2[j,i] = (z[j+1,i] - 2*z[j,i] + z[j-1,i])/dy**2

    zeta = d2zdx2 + d2zdy2

    return zeta

            
     
#
# define the function that returns the wind vector at given lat and lon and time step,
# as called by the euler() or huen() integration methods
def model_wind(lat, lon, time):

#
    # the interpolators are created in the main program, prior to calling this method
#    print("model_wind lat lon:", lat, lon)
    new_wind = wind_interpolator((lat,lon))
#    print("model_wind new wind:", new_wind)
    if new_wind is np.nan:
        return np.nan
    new_u = new_wind.real
    new_v = new_wind.imag
    return vel(new_u, new_v)  # this is the wind vector interpolated to the point
#

gravity = 9.81 # m/s/s
earth_radius = 6371.e3 # meters
# Create a grid of latitude and longitude values
# 1 degree grids
grid_spacing = 1. # 1 degree grid spacing
grid_spacing_rad = grid_spacing * np.pi / 180. # grid spacing in radians
min_lon = -160
max_lon = -20
min_lat = 10
max_lat = 70
nlat = max_lat - min_lat
#nlat, nlon = 140, 360
nlon = max_lon - min_lon
#lat_lin = np.linspace(-70, 70, nlat)
lat_lin = np.linspace(min_lat, max_lat, nlat)
mu = np.sin(lat_lin * np.pi/180.) # used in equation for geostrophic east-west wind
dmu = np.gradient(mu) 
#lon_lin = np.linspace(-180, 180, nlon)
lon_lin = np.linspace(min_lon, max_lon, nlon)
lon, lat = np.meshgrid(lon_lin, lat_lin)
phi_rad = np.array(lat_lin * np.pi/180.)  # latitude in radians
lambda_rad = np.array(lon_lin * np.pi/180.)  # longitude in radians

# Define the geopotential field with two ridges in the northern hemisphere and two ridges in the southern hemisphere
# north pacific
def prescribe_z(lon,lat):
    decay = 10
    r1 = np.exp(-((lat - 30) / decay) ** 2 - ((lon + 120) / decay) ** 2)

#    r2 = np.exp(-((lat + 30) / decay) ** 2 - ((lon + 120) / decay) ** 2)
    # north atlantic low
    r3 = - np.exp(-((lat - 30) / decay) ** 2 - ((lon + 80) / decay) ** 2)
#    r4 = np.exp(-((lat + 30) / decay) ** 2 - ((lon + 10) / decay) ** 2)


    return ((r1 + r3) * 100 + 5000) * gravity



def prescribe_winds():
# Compute the geostrophic winds
    dzdy = np.gradient(geopot, axis=0) / np.gradient(lat * np.pi / 180, axis=0)
    #
    #
    # calculate dzdy using constant difference argument of 1 deg resolution
    #
    dzdy = np.gradient(geopot,grid_spacing_rad, axis=0) 

    f = 2 * np.pi / 86400 * np.sin(lat * np.pi / 180)

    #
    #
    # dzdx is change of phi with respect to lambda
    #
    dzdx = np.gradient(geopot, grid_spacing_rad, axis=1)

    ug = - dzdy / (f* earth_radius)
    vg = dzdx/ (f * earth_radius * np.cos(lat * np.pi/180)) 

    wind_vector = np.vectorize(complex)(ug, vg)

    #
    return wind_vector

def plot_trajectories(init_lons, init_lats, deltat, nsteps):

    # plot trajectories
    fig4 = plt.figure(figsize=(12,8))


    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='black')
    m.drawcountries(linewidth=0.5, color='black')
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])

    hours = deltat/3600.  # number of hours in time step
    elapsed = nsteps * hours
    title_text = "Particle Trajectories, Time step = {:.1f} hours, No. steps: {:d}".format(hours, nsteps)
    plt.title(title_text)
# the colors of the trajectories cycle through the following list
    colors = ['black', 'red', 'blue', 'green','grey','orange', 'purple']
    for j in range(len(init_lons)):
        parcel_trajectory = huen(model_wind, init_lats[j], init_lons[j], deltat, nsteps)


        lons = parcel_trajectory.lons()
        lats = parcel_trajectory.lats()
#        print("Traj lats: ", lats)
#        print("Traj lons: ", lons)

#

        color_index = j % len(colors)
        line_color = colors[color_index]
        for i in range(len(lats)-1):
            dx = lons[i+1] - lons[i]
            dy = lats[i+1] - lats[i]
            plt.arrow(lons[i], lats[i], dx, dy, \
                      length_includes_head=True, head_length=0.8, head_width=0.8, color=line_color)

    plt.savefig('trajectories.png')
    plt.show()

#
# generate geopotential field on the lat lon grid
#
geopot = prescribe_z(lon,lat)

winds = prescribe_winds()
wsize =winds.shape

## Create a new figure
fig = plt.figure(figsize=(12, 6))

# Create a new map projection
m = Basemap(projection='cyl', llcrnrlat=min_lat, urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon)

# Draw the continents and coastlines
m.drawcoastlines(linewidth=0.5,color='white')
m.drawcountries(linewidth=0.5, color='white')

# Draw the geopotential field
x, y = m(lon, lat)
m.contourf(x, y, geopot, cmap='jet', levels=50)

# Draw the geostrophic wind vectors
ug = winds.real
vg = winds.imag


m.quiver(x[::5, ::5], y[::5, ::5], ug[::5, ::5], vg[::5, ::5], scale=1000, color='white')

m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
# Add a colorbar and title
plt.colorbar(label='Geopotential')
plt.title('Geopotential and Geostrophic Wind Vectors')
plt.savefig('winds.png')

plt.show()


#define methods to interpolate u and v wind components


wind_interpolator = scipy.interpolate.RegularGridInterpolator((lat_lin, lon_lin),\
                                                              winds, method='linear',bounds_error=False)


f = 2 * np.pi / 86400 * np.sin(lat * np.pi / 180)

# Create a separate plot for the geostrophic wind speed
fig2 = plt.figure(figsize=(12, 8))

# Draw the continents and coastlines in white
m.drawcoastlines(linewidth=0.5, color='white')
m.drawcountries(linewidth=0.5, color='white')

m.contourf(x,y,np.abs(winds), cmap='jet',levels=50)
m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])

# Draw the geostrophic# Show the plot
# Add a colorbar and title
plt.colorbar(label='Wind speed')
plt.title('Geostrophic Wind Speed')
plt.savefig("windspeed.png")
plt.show()

#
# plot the absolute voriticity
#


coslat = np.cos(lat * np.pi/180)
dvdx = np.gradient(vg, axis=0)/np.gradient(lon * np.pi/180, axis=1) / (earth_radius*coslat)
dudy = np.gradient(ug,axis=1)/np.gradient(lon * np.pi/180, axis=1) / (earth_radius*coslat) 

#
# use laplacian to get vorticity
#

dzdphi = np.gradient(geopot,grid_spacing_rad, axis=0)
dzdlambda = np.gradient(geopot,grid_spacing_rad, axis=1)
d2zdphi2 = np.gradient(dzdphi,grid_spacing_rad, axis=0)
d2zdlambda2 = np.gradient(dzdlambda, grid_spacing_rad, axis=1)

denom = 1./(earth_radius *np.cos(lat * np.pi/180))**2
relative_vorticity = d2zdphi2/earth_radius**2 + np.multiply(d2zdlambda2 ,denom)
absolute_vorticity = relative_vorticity + f

#
# recompute vorticity manually using centered differences
# results agree to within 2%.
#
#relative_vorticity = get_zeta(geopot,  lon, lat)
print("max rel vort: ", np.max(relative_vorticity))
print("min relativevort: ", np.min(relative_vorticity))
      
fig3 = plt.figure(figsize=(12,8))


# Draw the continents and coastlines in white                                                                            
m.drawcoastlines(linewidth=0.5, color='white')
m.drawcountries(linewidth=0.5, color='white')

m.contourf(x,y,relative_vorticity, cmap='jet',levels=50)
# Add a colorbar and title                                                                                      
m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
plt.colorbar(label='vorticity')
plt.title('relative vorticity')
plt.savefig("vorticity.png")
plt.show()
#
npart = 10 # number of parcel trajectories
nt = 10 # number of time steps
#

# Plotting the trajectories
deltat = 12 * 3600 # 12 hour time steps
nsteps = 10 # number of times to integrate over

#
# specify the initial lat lon points for trajectories
#
init_lats =  [40,40,45,40,40, 30, 30, 25, 30, 30, 45, 45, 25]
init_lons = [-60, -80,-100, -120,-140, -60, -80,-100, -120,-140, -120, -80, 80]
plot_trajectories(init_lons,init_lats, deltat, nsteps)
