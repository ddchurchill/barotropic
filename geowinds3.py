import cmath # used for complex numbers - storing wind vectors
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.interpolate
import xarray as xr

from velocity import Velocity as vel
from trajectory import Trajectory 
from traject import euler
from traject import huen
from diff import centered_diff
# constants

omega = 7.292e-5  # Earth's angular velocity in rad/s

def wind_and_vorticity(lon_grid, lat_grid, z_field, lat_spacing, lon_spacing):
    
    # Compute the numerical derivatives of the height field
    dz_dlat, dz_dlon = np.gradient(z_field, lat_spacing, lon_spacing)
    
    # Convert the derivatives to spherical coordinates
    dz_dphi = dz_dlat / (earth_radius * np.cos(np.deg2rad(lat_grid)))
    dz_dtheta = dz_dlon / earth_radius
    
    # Calculate the geostrophic wind components
    U_g = -(1 / f) * dz_dphi
    V_g = (1 / f) * dz_dtheta

    # Compute the numerical derivatives of the geostrophic wind components
    dUg_dlat, dUg_dlon = np.gradient(U_g, lat_spacing, lon_spacing)
    dVg_dlat, dVg_dlon = np.gradient(V_g, lat_spacing, lon_spacing)
    
    # Convert the derivatives to spherical coordinates
    dUg_dphi = dUg_dlat / (earth_radius * np.cos(np.deg2rad(lat_grid)))
    dVg_dtheta = dVg_dlon / earth_radius
    
    # Calculate the vorticity
    vorticity = dVg_dtheta - dUg_dphi

    print("numerical Max vorticity:", np.max(vorticity))
    return U_g, V_g, vorticity


def get_zeta1(geopot):
    #
    # use laplacian to get vorticity
    #

#    dzdphi = np.gradient(geopot,grid_spacing_rad, axis=0)
 #   dzdlambda = np.gradient(geopot,grid_spacing_rad, axis=1)
#    d2zdphi2 = np.gradient(dzdphi,grid_spacing_rad, axis=0)
# try passing lat and lon as 1-D arrays.  This works correctly
    # to within 1% of other methods. 
    dzdphi = np.gradient(geopot,lat_lin, axis=0)
    dzdlambda = np.gradient(geopot,lon_lin, axis=1)
    d2zdphi2 = np.gradient(dzdphi,lat_lin, axis=0)
    d2zdlambda2 = np.gradient(dzdlambda, lon_lin, axis=1)

    denom = 1./(earth_radius *np.cos(np.radians(lat)))**2
    relative_vorticity = d2zdphi2/earth_radius**2 + \
        np.multiply(d2zdlambda2 ,denom)

    print("get_zeta1, max vort:", np.max(relative_vorticity))
    return relative_vorticity

def get_zeta(z, x, y):
    """
    get_zeta - compute relative vorticity from geopotential height, z, with input longitudes, x ,

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
    print("get zeta, max vort:", np.max(zeta))
    return zeta

            
     
#
# define the function that returns the wind vector at given lat and lon and time step,
# as called by the euler() or huen() integration methods
# TODO: time is not used here yet.
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
nlat = max_lat - min_lat +1

nlon = max_lon - min_lon +1

lat_lin = np.linspace(min_lat, max_lat, nlat)
#print("lat_lin: ",lat_lin)
mu = np.sin(lat_lin * np.pi/180.) # used in equation for geostrophic east-west wind
dmu = np.gradient(mu) 

lon_lin = np.linspace(min_lon, max_lon, nlon)
lon, lat = np.meshgrid(lon_lin, lat_lin)
phi_rad = np.array(lat_lin * np.pi/180.)  # latitude in radians
lambda_rad = np.array(lon_lin * np.pi/180.)  # longitude in radians
lat_rad = np.radians(lat)  # convert degrees to radians
lon_rad = np.radians(lon)
f = 2 * omega * np.sin(lat_rad)  # Coriolis parameter


def algebraic_derivatives(lon, lat, decay=10, center_lat=40, gravity=9.81):
    r1 = np.exp(-((lat - center_lat) / decay) ** 2 - ((lon + 120) / decay) ** 2)
    r3 = -np.exp(-((lat - center_lat) / decay) ** 2 - ((lon + 80) / decay) ** 2)
 
    dr1_dlon = -200 * gravity * (lon + 120) / decay ** 2 * r1
    dr1_dlat = -200 * gravity * (lat - center_lat) / decay ** 2 * r1

    dr3_dlon = - 200* gravity * (lon + 80) / decay ** 2 * r3
    dr3_dlat = - 200* gravity * (lat - center_lat) / decay ** 2 * r3

    dz_dlon = (dr1_dlon + dr3_dlon) 
    dz_dlat = (dr1_dlat + dr3_dlat) 


# code for this method below is from chatgpt
# Assuming lat, lon, center_lat, decay, and gravity are defined
    lat_rad = np.radians(lat)
    center_lat_rad = np.radians(center_lat)
    lon_rad = np.radians(lon)

    
    R = 6371e3  # Earth's radius in meters

    r1 = np.exp(-((lat_rad - center_lat_rad) / decay) ** 2 - ((lon_rad + np.radians(120)) / decay) ** 2)
    r3 = -np.exp(-((lat_rad - center_lat_rad) / decay) ** 2 - ((lon_rad + np.radians(80)) / decay) ** 2)

# Partial derivative of z w.r.t lat (in radians)
    dz_dlat_rad = (2 * (lat_rad - center_lat_rad) / decay**2 * (r1 - r3)) * 100 * gravity

# Convert the derivative to degrees
    dz_dlat = dz_dlat_rad # * (180 / np.pi)

# Partial derivative of z w.r.t lon (in radians)
    dz_dlon_rad = (2 * (lon_rad + np.radians(120)) / decay**2 * r1 - 2 * (lon_rad + np.radians(80)) / decay**2 * r3) * 100 * gravity

# Convert the derivative to degrees
    dz_dlon = dz_dlon_rad # * (180 / np.pi)

# Compute the distance between points along the latitude and longitude lines
    dlat_dist = R * np.radians(dz_dlat)
    dlon_dist = R * np.cos(lat_rad) * np.radians(dz_dlon)

    return dz_dlon, dz_dlat

#print(f"dz/dlon: {dz_dlon}, dz/dlat: {dz_dlat}")

# Define the geopotential field with two ridges in the northern hemisphere and two ridges in the southern hemisphere
# north pacific
def ridge_and_trough(lon,lat):
    decay = 10
    center_lat = 40
    r1 = np.exp(-((lat - center_lat) / decay) ** 2 - ((lon + 120) / decay) ** 2)


    r3 = - np.exp(-((lat - center_lat) / decay) ** 2 - ((lon + 80) / decay) ** 2)

    z = ((r1 + r3) * 100 + 5000) * gravity

    return ((r1 + r3) * 100 + 5000) * gravity

def zonal_wind(lon,lat):  # make an east-west wind

     # make a rsmp from 45 to 35 degrees

# Define the number of longitude and latitude points
#    n_lon, n_lat = max_lon - min_lon , max_lat - min_lat
#    print("nlat nlon:", n_lat, n_lon)
# Create a 2D array of zeros with the given dimensions
    z = np.zeros((nlat, nlon))
    deltaz = 105.43 # meters height difference
    lowz = 5000 # starting height in meters
# Fill the height_data array with the desired height values
    for i in range(nlat):
        lat = max_lat - i
        if lat >= 45:
            z[nlat - i -1, :] = lowz
        elif 35 <= lat < 45:
            z[nlat -i -1, :] = lowz + deltaz * (45 - lat) / (45 - 35)
        else:
            z[nlat -i -1, :] = lowz + deltaz
    return z

def north_wind(lon, lat):
# Create a 2D array of zeros with the given dimensions
    z = np.zeros((nlat, nlon))
    max_z = 5600
    min_z = 5000
    left_lon = -120
    right_lon = -100
# Fill the height_data array with the desired height values
    for i in range(nlon):
        long = lon[0,i]
        if long < -120:
            z[:, i] = min_z
        elif -120 <= long < -100:
            z[:, i] = min_z + (max_z - min_z) * (long - left_lon) / (right_lon - left_lon)
        else:
            z[:, i] = max_z

    return z

# get barotropic forecast data from file
def baro_fcst(forecast_init_time):
    step = fc_ds_baro.step[0]
    print("forecast time:", forecast_init_time)
    baro = fc_ds_baro.sel({
        "time": forecast_init_time,
        "step": step,
        "lat": slice(max_lat, min_lat),
        "lon": slice(180+min_lon, 180+max_lon)
    })['z500'].values
    
    print("baro min height:", np.min(baro))
    print("baro max height:", np.max(baro))      

    # check the lat and lon range
    lats = fc_ds_baro.sel({
        "time": forecast_init_time,
        "step": step,
        "lat": slice(max_lat, min_lat),
        "lon": slice(180+min_lon, 180+max_lon)
    })['lat'].values
    lons = fc_ds_baro.sel({
        "time": forecast_init_time,
        "step": step,
        "lat": slice(max_lat, min_lat),
        "lon": slice(180+min_lon, 180+max_lon)
    })['lon'].values
#    print("Baro_fcst lats:", np.flip(lats))
#    print("Baro_fcst lons:", np.flip(lons))
    #flip the data around the latitude axis, as it is stored in the
    # file upside down compared to my usage.
    baro = np.flip(baro,axis=0)
    return baro

def prescribe_winds1():
#
    # use the algrebaric derivative of the z field to get
    # the winds
    dzdx, dzdy = algebraic_derivatives(lon, lat)

    ug = - dzdy / f /earth_radius
    vg = dzdx/ f /earth_radius
    wind_vector = np.vectorize(complex)(ug, vg)
    speed = np.sqrt(ug*ug + vg*vg)
    print("prescribe winds1: max wind speed:", np.max(speed))
    x = np.cos(lat_rad)*earth_radius*np.cos(lon_rad)
    y = earth_radius * lat_rad
    dvdx, dydy = centered_diff(vg, x, y)
    dudx, dudy = centered_diff(ug, x, y)
    zeta = dvdx - dudy
    print("prescribe winds1: max vort: ",np.max(zeta))
    return wind_vector, zeta

def prescribe_winds2():
    # use centered differences routine to get winds
    print("prescribe_winds2:")
    x = np.cos(lat_rad)*earth_radius*lon_rad
    y = earth_radius * lat_rad
    dzdx , dzdy =  centered_diff(geopot, x, y)
    #
    ug = - gravity *dzdy/f
    vg = gravity * dzdx /f
    wind_vector = np.vectorize(complex)(ug,vg)
    speed = np.sqrt(ug*ug + vg*vg)
    dvdx, dvdy = centered_diff(vg, x, y)
    dudx, dudy = centered_diff(ug, x, y)
    zeta = dvdx - dudy
    return wind_vector, zeta, speed

def prescribe_winds():
# Compute the geostrophic winds
    x = np.cos(lat_rad)*earth_radius*lon_rad
    y = earth_radius * lat_rad
    dzdx = np.gradient(geopot, axis=1)/np.gradient(x, axis=1)
    dzdy = np.gradient(geopot, axis=0)/np.gradient(y,axis=0)
    ug = - gravity *dzdy/f
    vg = gravity * dzdx /f
    wind_vector = np.vectorize(complex)(ug,vg)
    speed = np.sqrt(ug*ug + vg*vg)
    dvdx = np.gradient(vg, axis=1)/np.gradient(x,axis=1)
    dudy = np.gradient(ug, axis=0)/np.gradient(y,axis=0)
    zeta = dvdx - dudy
    return wind_vector, zeta, speed
def plot_winds(wind_data):
        ## Create a new figure
    fig = plt.figure(figsize=(12, 8))
    

    # Draw the continents and coastlines
    m.drawcoastlines(linewidth=0.5,color='white')
    m.drawcountries(linewidth=0.5, color='white')

    # Draw the geopotential field
    x, y = m(lon, lat)
    m.contourf(x, y, geopot, cmap='jet', levels=50)
    
    # Draw the geostrophic wind vectors
    ug = wind_data.real
    vg = wind_data.imag
    
    
    m.quiver(x[::5, ::5], y[::5, ::5], ug[::5, ::5], vg[::5, ::5], \
         scale=2000, color='white')
    
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    # Add a colorbar and title
    plt.colorbar(label='Geopotential')
    plt.title('Geopotential and Geostrophic Wind Vectors, ' + dt_str)
    plt.savefig('winds'+dt_str+'.png')

    plt.show()

def plot_trajectories(deltat, nsteps):

    # plot trajectories
    fig4 = plt.figure(figsize=(12,8))


    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='black')
    m.drawcountries(linewidth=0.5, color='black')
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])

    hours = deltat/3600.  # number of hours in time step
    elapsed = nsteps * hours
    title_text = "Trajectories, Time step = {:.1f} hours, No. steps: {:d}. ".format(hours, nsteps)
    plt.title(title_text + dt_str)
# the colors of the trajectories cycle through the following list
    colors = ['black', 'red', 'blue', 'green','grey','orange', 'purple']
#
# compute trajectories every 10 degrees in the domain
    lat_range = range(min_lat, max_lat +1, 10)
    lon_range = range(min_lon, max_lon +1, 10)
    n_trajectories = len(lat_range) * len(lon_range)
    color_index = 0
#
# repeat for each latitude and longitude in range
#
# create an empty list of trajectories
    trajectory_list = []

    index = 0
    for lat0 in lat_range:
        for lon0 in lon_range:
            # compute the trajectory from the model winds
            trajectory = \
                huen(model_wind, lat0, lon0, deltat, nsteps)
            # append the parcel to the trajectory list
            trajectory_list.append(trajectory)

            # here is where we compute the vorticity along the trajectory
            # iterate through the trajectory
            # TODO
            lons = trajectory.lons()
            lats = trajectory.lats()
                
            # plot the trajectory
            color_index = color_index + 1
            line_color = colors[color_index % len(colors)]
            for i in range(trajectory.length):
                dx = lons[i+1] - lons[i]
                dy = lats[i+1] - lats[i]
                plt.arrow(lons[i], lats[i], dx, dy, \
                      length_includes_head=True, head_length=0.8, head_width=0.8, color=line_color)
                # determine the vorticity at this lat/lon point
                vort = vort_interpolator((lats[i],lons[i]))
                # append it to the vorticity array in the trajectory
                trajectory.vort.append(vort)
#
# TODO: Analytics on changes in vorticity
#

    plt.savefig('trajectories'+dt_str+'.png')
    plt.show()
#

def plot_speed(wind_data):
        # Create a separate plot for the geostrophic wind speed
    fig2 = plt.figure(figsize=(12, 8))
    
    # Draw the continents and coastlines in white
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')
    
    x, y = m(lon, lat)

    m.contourf(x,y,np.abs(wind_data), cmap='jet',levels=50)
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    
    # Draw the geostrophic# Show the plot
    # Add a colorbar and title
    plt.colorbar(label='Wind speed')
    plt.title('Geostrophic Wind Speed ' + dt_str)
    plt.savefig('windspeed' + dt_str + '.png')
    plt.show()

def plot_vort(vort):
        
#    coslat = np.cos(lat * np.pi/180)

    absolute_vort = vort + f

    #
      
    fig3 = plt.figure(figsize=(12,8))
    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')

    x, y = m(lon, lat)

    m.contourf(x,y,absolute_vort, cmap='jet',levels=50)
    # Add a colorbar and title                                                                                      
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
    plt.colorbar(label='vorticity')
    plt.title('absolute vorticity ' + dt_str)
    plt.savefig("abs_vort"+dt_str+".png")
    plt.show()


# main program:
# Open the NetCDF file
fc_ds_baro = xr.open_dataset("forecasts_baro-spharm_2007-04-15_2007-04-18.nc")
fc_ds_baro = fc_ds_baro.roll({"lon": 180})
# Create a new map projection
m = Basemap(projection='cyl', llcrnrlat=min_lat, urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon)

for time_stamp in fc_ds_baro['time'].values:

          

#dt_str = np.datetime64.strftime(timestamp, '%Y-%m-%d %H:%M:%S')
# set unit=D for days, =s for seconds
    dt_str = np.datetime_as_string(time_stamp, unit='s')
    print("time stamp ", dt_str)
    geopot = baro_fcst(time_stamp)  # read in height field from data file.

    # generate geopotential field on the at lon grid
    #
    #geopot = ridge_and_trough(lon,lat)
    #geopot = zonal_wind(lon,lat)  # create a westerly wind only
    #geopot = north_wind(lon,lat)  # create a westerly wind only
    # Calculate the geostrophic wind components and vorticity
    #U_g, V_g, vorticity = wind_and_vorticity(lon_grid, lat_grid, \
    # dz_dphi_algebraic, dz_dtheta_algebraic, lat_spacing, lon_spacing)

    winds2, zeta2, speed2 = prescribe_winds2() # uses my differences code
    #
    # use the gradient calls
    #
    winds, zeta, speed = prescribe_winds() # numpy gradient method
    relative_vorticity = zeta2 # use my differences code
    rms_speed = np.sqrt(np.mean(np.square(speed2)))
    rms_speed_diff = np.sqrt(np.mean(np.square(speed2 - speed)))
    print("max wind speed centered:", np.max(speed2))
    print("Max wind speed gradient:", np.max(speed))
    print("rms speed: ", rms_speed)
    print("rms speed diff:", rms_speed_diff)
    print("max vort: ",np.max(relative_vorticity))
    print("min vort: ",np.min(relative_vorticity))
    rms_diff = np.sqrt( np.mean(np.square(zeta - zeta2)))
    rms_vort = np.sqrt( np.mean(np.square(zeta)))
    print("vort rms diff: ", rms_diff)
    print("vort rms: ", rms_vort)
    print("vort: rel diff %:", rms_diff/rms_vort*100)
    #
    # print wind speed and vorticity by latitude
    #
    #    print("lat, speed, vorticity")
    #    for j in range(1,59):
    #        print(lat[j,0], speed[j,0], zeta[j,0])

    plot_winds(winds2) # winds2 is manual centered differences


    #define methods to interpolate u and v wind components
    # this gets called by model_winds() which is called by
    # plot_trajectories()
    wind_interpolator = scipy.interpolate.RegularGridInterpolator((lat_lin, lon_lin),\
            winds2, method='linear',bounds_error=False)
    # create an interpolator for vorticity
    absolute_vort = zeta2 + f
    vort_interpolator = scipy.interpolate.RegularGridInterpolator((lat_lin, lon_lin),\
            absolute_vort, method='linear',bounds_error=False)

    plot_speed(winds2)
    

    #
    # plot the relative voriticity
    #
    plot_vort(zeta2)    


    # Plotting the trajectories
    deltat = 1 * 3600 # 1 hour time steps
    nsteps = 12 # number of times to integrate over
    #
    plot_trajectories(deltat, nsteps)
