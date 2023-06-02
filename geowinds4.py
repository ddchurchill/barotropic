from scipy.misc import derivative
import copy
import matplotlib.gridspec as gridspec
import os
import cmath # used for complex numbers - storing wind vectors
import math
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from mpl_toolkits.basemap import Basemap
from trajectory_v2 import Trajectory_v2
from huen_v4 import *
import scipy.interpolate
import xarray as xr

from velocity import Velocity as vel
#from trajectory import Trajectory 
#from traject import euler
#from traject import huen
from diff import *
from traject_constants import * 
from plot_traject_arrows import *
# constants

def get_timestamp(start_time, time_step):
    time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_step))
    dt_str = np.datetime_as_string(time_stamp, unit='s')
    return dt_str, time_stamp


def wind_and_vorticity(lon_grid, lat_grid, z_field, lat_spacing, lon_spacing):
    
    # Compute the numerical derivatives of the height field
    dz_dlat, dz_dlon = np.gradient(z_field, lat_spacing, lon_spacing)
    
    # Convert the derivatives to spherical coordinates
    dz_dphi = dz_dlat / (EARTH_RADIUS * np.cos(np.deg2rad(lat_grid)))
    dz_dtheta = dz_dlon / EARTH_RADIUS
    
    # Calculate the geostrophic wind components
    U_g = -(1 / f) * dz_dphi
    V_g = (1 / f) * dz_dtheta

    # Compute the numerical derivatives of the geostrophic wind components
    dUg_dlat, dUg_dlon = np.gradient(U_g, lat_spacing, lon_spacing)
    dVg_dlat, dVg_dlon = np.gradient(V_g, lat_spacing, lon_spacing)
    
    # Convert the derivatives to spherical coordinates
    dUg_dphi = dUg_dlat / (EARTH_RADIUS * np.cos(np.deg2rad(lat_grid)))
    dVg_dtheta = dVg_dlon / EARTH_RADIUS
    
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

    denom = 1./(EARTH_RADIUS *np.cos(np.radians(lat)))**2
    relative_vorticity = d2zdphi2/EARTH_RADIUS**2 + \
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
            dx = EARTH_RADIUS * np.cos(np.radians(y[j,i])) * np.radians( dlon)
            d2zdx2[j,i] = (z[j,i+1] - 2*z[j,i] + z[j,i-1])/dx**2
 
            dlat = (y[j+1,i] - y[j-1,i])/2.
#            print("dlat:", dlat)
            dy = EARTH_RADIUS * np.radians(dlat)
            d2zdy2[j,i] = (z[j+1,i] - 2*z[j,i] + z[j-1,i])/dy**2

    zeta = d2zdx2 + d2zdy2
    print("get zeta, max vort:", np.max(zeta))
    return zeta

            
     
#
# interp_data - interpolate data from dataset
# input: field - name of field to get
# lat0: latitude of point
# lon0: longitude of point
# time0: time stamp of point
#
def interp_data(field, lat0, lon0, time0):

    interpolated_value = fc_ds_baro[field].interp(
    lon=lon0,
    lat=lat0,
    time=time0,
    method='linear')
    return interpolated_value
#
# interp_wind - return wind interpolated to given lat, lon and timestamp
#
def interp_wind(lat0, lon0, time0):
    value = fc_ds_baro['winds'].interp( \
         lon=lon0, \
         lat=lat0, \
         time=time0, \
         method='linear')
    return value


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

lon_lin = np.linspace(min_lon, max_lon, nlon)
lon, lat = np.meshgrid(lon_lin, lat_lin)
phi_rad = np.array(lat_lin * np.pi/180.)  # latitude in radians
lambda_rad = np.array(lon_lin * np.pi/180.)  # longitude in radians
lat_rad = np.radians(lat)  # convert degrees to radians
lon_rad = np.radians(lon)
f = 2 * OMEGA * np.sin(lat_rad)  # Coriolis parameter


# Define the geopotential field with two ridges in the northern hemisphere and two ridges in the southern hemisphere
# north pacific
def ridge_and_trough():
    decay = 10
    center_lat = 40
    r1 = np.exp(-((lat - center_lat) / decay) ** 2 - ((lon + 120) / decay) ** 2)


    r3 = - np.exp(-((lat - center_lat) / decay) ** 2 - ((lon + 80) / decay) ** 2)

    z = ((r1 + r3) * 1000 + 5000) 

    return z

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
#
# overwrite the dataset height values with test data
#
def south_wind_v2(z, lats, lons):
    max_z = 5600
    min_z = 5000
    left_lon = -120
    right_lon = -100
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            if lon < -120:
                value = min_z
            elif -120 <= lon < -100:
                value = min_z + (max_z - min_z) * \
                    (lon - left_lon) / (right_lon - left_lon)
            else:
                value = max_z
            z[i][j] = value
    return z
def north_wind_v2(z, lats, lons):
    max_z = 5600
    min_z = 5000
    left_lon = -120
    right_lon = -100
    scale = - (max_z - min_z)/(right_lon - left_lon)
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            if lon < left_lon:
                z[i,j] = max_z
            elif left_lon <= lon < right_lon:
                z[i,j] = min_z + scale * (lon - right_lon)
            else:
                z[i,j] = min_z
    return z


    
# get barotropic forecast data from file
def baro_fcst(forecast_init_time, step):
    print("forecast time:", forecast_init_time)
    baro = fc_ds_baro.sel({
        "time": forecast_init_time,
        "step": step,
        "lat": slice(max_lat, min_lat),
        "lon": slice(180+min_lon, 180+max_lon)
    })['z500'].values
    
#    print("baro min height:", np.min(baro))
#    print("baro max height:", np.max(baro))      

    # check the lat and lon range
#    lats = fc_ds_baro.sel({
#        "time": forecast_init_time,
#        "step": step,
#        "lat": slice(max_lat, min_lat),
#        "lon": slice(180+min_lon, 180+max_lon)
#    })['lat'].values
#    lons = fc_ds_baro.sel({
#        "time": forecast_init_time,
#        "step": step,
#        "lat": slice(max_lat, min_lat),
#        "lon": slice(180+min_lon, 180+max_lon)
#    })['lon'].values
#    print("Baro_fcst lats:", np.flip(lats))
#    print("Baro_fcst lons:", np.flip(lons))
    #flip the data around the latitude axis, as it is stored in the
    # file upside down compared to my usage.
    baro = np.flip(baro,axis=0)
#
    return baro
def prescribe_winds3(): 
    """
    prescriber_winds3 ;
    use mu = sin(theta) as north-south cooridinate
    use gradient calls to diferentiate
    """
    mu = np.sin(lat_rad) # 2-D array north-south component
    mu1 = np.sin(np.deg2rad(lat_lin)) # 1-D array

    phi = GRAVITY * geopot
    dphidmu = np.gradient(phi, mu1, axis=0)
    dphidmu1 = (1 - mu**2)*dphidmu
    d2phidmu2 = np.gradient(dphidmu1, mu1,  axis=0)

# east-west component 
#
    lambdax = np.cos(lat_rad)* EARTH_RADIUS * lon_rad
    dlambdax = np.gradient(lambdax, axis=1)
# first derivative w/r longitude
    dphidlambda = np.gradient(phi,axis=1)/dlambdax
# 2nd derivative w/r longitude
    d2phidlambda2 = np.gradient(dphidlambda,axis=1)/dlambdax 
# scaled by 1/cosine squared
    d2phi1 = d2phidlambda2 / ( 1 - mu**2) 
#
# Add north-south and east-west second derivatives
# Note need to divide by f, otherwise too small by 10^-1!  5/9/23
    zeta =  (d2phidmu2 + d2phi1) / EARTH_RADIUS**2 /f
#
# compute geostrophic winds
#
    ug = -  dphidmu * np.cos(lat_rad)/EARTH_RADIUS /f
    vg =  dphidlambda /f

    # compute wind speed
    speed = np.sqrt(ug*ug + vg*vg)
    return ug, vg, zeta,speed

def prescribe_winds2():
    """
    prescribe_winds2 - uses latitude as north-south cooridnate
    uses manual centered differences to differentiate
    """
    x = np.cos(lat_rad)*EARTH_RADIUS*lon_rad
    y = EARTH_RADIUS * lat_rad
    dzdx , dzdy =  centered_diff(geopot, x, y)
    #
    ug = - GRAVITY *dzdy/f
    vg = GRAVITY * dzdx /f
#    wind_vector = np.vectorize(complex)(ug,vg)
    speed = np.sqrt(ug*ug + vg*vg)
    dvdx, dvdy = centered_diff(vg, x, y)
    dudx, dudy = centered_diff(ug, x, y)
    zeta = dvdx - dudy
    return ug, vg, zeta, speed

def prescribe_winds(): 
    """
    prescribe_winds()
    uses latitude as north-south coordinate.
    uses gradient calls for derivatives˜
    """

    x = np.cos(lat_rad)*EARTH_RADIUS*lon_rad
    y = EARTH_RADIUS * lat_rad
    dzdx = np.gradient(geopot, axis=1)/np.gradient(x, axis=1)
    dzdy = np.gradient(geopot, axis=0)/np.gradient(y,axis=0)
    ug = - GRAVITY *dzdy/f
    vg = GRAVITY * dzdx /f
    wind_vector = np.vectorize(complex)(ug,vg)
    speed = np.sqrt(ug*ug + vg*vg)
    dvdx = np.gradient(vg, axis=1)/np.gradient(x,axis=1)
    dudy = np.gradient(ug, axis=0)/np.gradient(y,axis=0)
    zeta = dvdx - dudy
#
    return wind_vector, zeta, speed
def plot_winds(wind_data):
        ## Create a new figure
    fig = plt.figure(figsize=(12, 8))
    

    # Draw the continents and coastlines
    m.drawcoastlines(linewidth=0.5,color='white')
    m.drawcountries(linewidth=0.5, color='white')

    # Draw the geopotential field
    x, y = m(lon, lat)
    m.contourf(x, y, geopot, cmap='jet', levels=30, vmin=5000., vmax=6000.)
    
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

def plot_winds_v2(dataset, time_index, time_str, showplot):

#    time_stamp = dataset['time'][time_index].values
#    steps = dataset['step'].values
#    time_step = steps[time_index]
#    print("time index is ", time_index)
#    start_time = dataset['time'].values
#    print("plotwinds2: start time:", start_time)
#    time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_index))
#    dt_str = np.datetime_as_string(time_stamp, unit='s')
#

    print("plot winds v2: time stamp ", time_str)
    steps = dataset['step'].values
    
#    print("plot winds v23: time = ", time_stamp)
#    time_str = np.datetime_as_string(time_stamp, unit='s')
    print("plot winds v2: time string:", time_str)
        ## Create a new figure
    fig = plt.figure(figsize=(12, 8))
    

    # Draw the continents and coastlines
    m.drawcoastlines(linewidth=0.5,color='white')
    m.drawcountries(linewidth=0.5, color='white')

    # Draw the geopotential field
    x, y = m(lon, lat)
    heights = dataset['z500'][time_index].values
    m.contourf(x, y, heights, cmap='jet', levels=30, vmin=5000., vmax=6000.)
    
    # Draw the geostrophic wind vectors
    ug = dataset['wind_u'][time_index].values
    vg = dataset['wind_v'][time_index].values
    
    
    m.quiver(x[::5, ::5], y[::5, ::5], ug[::5, ::5], vg[::5, ::5], \
         scale=2000, color='white')
    
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    # Add a colorbar and title
    plt.colorbar(label='Geopotential')
    plt.title('Geopotential and Geostrophic Wind Vectors, ' + time_str)
    plt.savefig('winds'+time_str+'.png')

    if showplot:
        plt.show()
    plt.close()
    
    
#
def compute_trajectories(dataset, start_time, deltat, nsteps):
    trajectory_file = 'trajectories.npz' # numpy compressed file
    if os.path.exists(trajectory_file):
        npdata = np.load(trajectory_file,allow_pickle=True)
        trajectory_list = npdata['trajectories']

        print("Read in trajectory data from ", trajectory_file)
    else: # compute the trajectories
        lat_range = range(min_lat +5, max_lat -4, 10)
        lon_range = range(min_lon + 5, max_lon -5, 10)
        hours = deltat/3600.  # number of hours in time step
        subtitle_text = "Time step = " + str(hours) + " hours"
#        elapsed = nsteps * hours
        trajectory_list = []
        index = 0
        start_step = 0 # starting time step for trajectory
        for lat0 in lat_range:
            for lon0 in lon_range:
                # compute the trajectory from the model winds
                trajectory = \
                    huen_v4(dataset, lat0, lon0, start_step, nsteps, deltat)
                # append the parcel to the trajectory list
                trajectory_list.append(trajectory)

        # save the data to a file.
        np.savez(trajectory_file, trajectories=trajectory_list)
                 
        print("Wrote trajectories to ", trajectory_file)
    return trajectory_list
#
def plot_one_trajectory(traj, filename, start_time, stop_time):
    """
    input: traj - a trajectory
           filename - name of file to write plot to
    """
    fig = plt.figure(figsize=(8, 10))
    ax = fig.add_subplot(211, projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(212)
    ax.set_extent([min_lon, max_lon, min_lat, max_lat])
    # Draw the continents and coastlines in white
    ax.coastlines(linewidth=0.5, color='black')
    ax.add_feature(cartopy.feature.BORDERS, linewidth=0.5, edgecolor='black')

    # Draw parallels and meridians
    parallels = range(min_lat, max_lat, 10)
    meridians = range(min_lon, max_lon, 10)
    ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5,
                 linestyle='--')
    ax.set_xticks(meridians, crs=ccrs.PlateCarree())
    ax.set_yticks(parallels, crs=ccrs.PlateCarree())
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    #
    # iterate throught the points in a trajectory,
    # plot an arrow head for each
    #
    colors = ['black', 'red', 'blue', 'orange','grey','green', 'purple']
#
    color_index = 0
    color_interval = 7
    for i, point in enumerate(traj.points):
#        print(i, point.lat, point.lon, point.timestamp)
#        i += 1
        # make the color of the trajecotry change every 6 time steps
        # and cycle round after running through the 7 colors
        color_index = int((i / color_interval)) % (color_interval + 1)  
        ax.arrow(point.lon, point.lat, point.dx, point.dy, \
		  length_includes_head=True, head_length=1.0, \
		  head_width=1.0, color=colors[color_index])

    # plot a timeline of vorticity values     

    ax2.set_ylim(0, 2.e-4)
    ax2.set_xlabel('Date Time')
    ax2.set_ylabel('Abs. Vorticity')
    ax2.set_xlim(start_time, stop_time)
    t = traj
    title = "start:  {:.2f},{:.2f}, end: {:.2f},{:.2f}".\
        format(t.start_lat, t.start_lon, t.last_lat, t.last_lon)
    start_str = np.datetime_as_string(start_time, unit='s')
    # the string for the last time needs to show the final time in
    # the trajectory, which may be shorter than nsteps, in those cases
    # where the trajectory goes on the edge of the plot.
    #
    stop_str =  np.datetime_as_string(t.stop_time, unit='s')
    title = title + "\n" + start_str + " to " + stop_str
    line1 = "Begin: {:.2f}, {:.2f}, ".format(t.start_lat,t.start_lon)
    line2 = "End: {:.2f}, {:.2f}, ".format(t.last_lat,t.last_lon)
    title = line1 + start_str + "\n" + line2 + stop_str
#    print(title)
    plt.title(title)

    subtitle = "Time step = {:.1f} hours".format(t.time_step/3600)
    plt.suptitle(subtitle)
        
#
# plot the timeline of vorticity
#
    ax2.grid(axis='y', linestyle='--')
    ax2.grid(axis='x', linestyle='--')

    time_steps = [p.timestamp for p in traj.points] # time steps in hours
    # convert to date/time stamps,
    time_list = []
    for i, step in enumerate(time_steps):
       time_str, timestamp = get_timestamp(start_time, step)
       time_list.append(timestamp)
       
    # get the list of vorticity values from this trajectory
    vort_list = [p.vort for p in traj.points]
    #   plot the time line versus the vorticity
#    ax2.plot(time_list, vort_list,color='red')
    color_index = 0
    for i in range(0, len(vort_list), color_interval):
        plt.plot(time_list[i:i+color_interval+1], \
                 vort_list[i:i+color_interval+1], \
                 color=colors[color_index])
        color_index += 1
        color_index %= color_interval
    plt.tight_layout() # this is needed to prevent overlapping figures.

    plt.savefig(filename)
    plt.close()
#    plt.show()
#
# save the plot for a file
#

    
def plot_trajectories(trajectories, start_time, last_time):
    """
    input:
        trajectory list
        start_time - timestamp of start of period
        last_time - timestamp of end of period
    """

    fig4 = plt.figure(figsize=(12,8))


    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='black')
    m.drawcountries(linewidth=0.5, color='black')
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])

# the colors of the trajectories cycle through the following list
    colors = ['black', 'red', 'blue', 'green','grey','orange', 'purple']
#
    color_index = 0
#    start_str = str(start_time)
#    last_str = str(last_time)
    start_str = np.datetime_as_string(start_time, unit='s')
    last_str = np.datetime_as_string(last_time, unit='s')

    title_text = "Trajectories, " + start_str + " to " + last_str
    print(title_text)
    plt.title(title_text)
#    plt.suptitle(subtitle_text)        
# plot the trajectories
#
    plot_traject_arrows(trajectories)

    #
    file_name = "traj_" + start_str + "-" + last_str + ".png"
    plt.savefig(file_name)
    plt.show()
    plt.close()

def plot_speed_v2(dataset, time_index, time_str, showplot):

    fig2 = plt.figure(figsize=(12, 8))
     
    # Draw the continents and coastlines in white
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')
    
    x, y = m(lon, lat)

    speed = dataset['speed'][time_index].values
    m.contourf(x,y,speed, cmap='jet',levels=30)
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    
    # Draw the geostrophic# Show the plot
    # Add a colorbar and title
    plt.colorbar(label='Wind speed')
#    time_stamp = dataset['time'][time_index].values
#    time_str = np.datetime_as_string(time_stamp, unit='s')
    plt.title('Geostrophic Wind Speed ' + time_str)
    plt.savefig('windspeed' + time_str + '.png')
    if showplot:
        plt.show()
    plt.close()
    

def plot_vort(vort):
        
    absolute_vort = vort + f
    print("plot vort: min abs vort:", np.min(absolute_vort))
    print("plot vort: max abs vort:", np.max(absolute_vort))
    #
      
    fig3 = plt.figure(figsize=(12,8))
    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')

    x, y = m(lon, lat)

#    m.contourf(x,y,absolute_vort, cmap='jet',levels=30,vmin=0.,vmax=1.e-3)
#    absolute_vort[absolute_vort >0.] = 0 # filter out positive values
    m.contourf(x,y,absolute_vort, cmap='jet',levels=30)
    # Add a colorbar and title                                                                                      
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
    plt.colorbar(label='vorticity')
    time_stamp = dataset['time'][time_index].values
    time_str = np.datetime_as_string(time_stamp, unit='s')
    plt.title('absolute vorticity ' + time_str)
    plt.savefig("abs_vort"+time_str+".png")
    if showplot:
        plt.show()
    plt.close()
    
def plot_vort_v2(dataset, time_index, time_str, showplot):

#    time_stamp = dataset['time'][time_index].values
#    time_str = np.datetime_as_string(time_stamp, unit='s')
    vort = dataset['abs_vorticity'][time_index].values
    print("plot vort: min abs vort:", np.min(vort))
    print("plot vort: max abs vort:", np.max(vort))

    #
      
    fig3 = plt.figure(figsize=(12,8))
    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')

    x, y = m(lon, lat)

#    m.contourf(x,y,absolute_vort, cmap='jet',levels=30,vmin=0.,vmax=1.e-3)
#    absolute_vort[absolute_vort >0.] = 0 # filter out positive values
    m.contourf(x,y,vort, cmap='jet',levels=16, vmin=0., vmax=0.0003)
    # Add a colorbar and title                                                                                      
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
    plt.colorbar(label='vorticity')
    plt.title('absolute vorticity ' + time_str)
    plt.savefig("abs_vort"+time_str+".png")
    if showplot:
        plt.show()
    plt.close()

def plot_rel_vort(vort): # plot relative voriticity
        
    fig3a = plt.figure(figsize=(12,8))
    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')

    x, y = m(lon, lat)

#    m.contourf(x,y,absolute_vort, cmap='jet',levels=30,vmin=0.,vmax=1.e-3)
# plot only the negative values of the relative vorticity to see where
# they are. 
#    vort[vort > .0] = 0
    m.contourf(x,y,vort, cmap='jet',levels=15)
    # Add a colorbar and title                                                                                      
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
    plt.colorbar(label='vorticity')
    plt.title('relative vorticity ' )
    plt.savefig("rel_vort"+dt_str+".png")
    plt.show()

def plot_rel_vort_v2(dataset, time_index, time_str, showplot):

        
#    time_stamp = dataset['time'][time_index].values
#    time_str = np.datetime_as_string(time_stamp, unit='s')
    vort = dataset['rel_vorticity'][time_index].values

    fig3a = plt.figure(figsize=(12,8))
    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='white')
    m.drawcountries(linewidth=0.5, color='white')

    x, y = m(lon, lat)

    m.contourf(x,y,vort, cmap='jet',levels=30)
    # Add a colorbar and title                                                                                      
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
    plt.colorbar(label='vorticity')
    plt.title('relative vorticity ' + time_str)
    plt.savefig("rel_vort"+time_str+".png")
    if showplot:
        plt.show()
    plt.close()
    
def plot_vort_advection(dataset, time_index, time_str, showplot):

        
    zeta = dataset['abs_vorticity'][time_index].values
    u = dataset['wind_u'][time_index].values
    v = dataset['wind_v'][time_index].values
    x = np.cos(lat_rad)*EARTH_RADIUS*lon_rad
    y = EARTH_RADIUS * lat_rad
    dzetadx = np.gradient(zeta, axis=1)/np.gradient(x, axis=1)
    dzetady = np.gradient(zeta, axis=0)/np.gradient(y,axis=0)

    vadv = -(u * dzetadx + v*dzetady)
    vadv[np.abs(vadv) < 1.e-16 ] = 0 # filter noise
    print("max value of vort advection is ", np.max(np.abs(vadv)))
    fig = plt.figure(figsize=(12,8))
    # Draw the continents and coastlines in white                                                                            
    m.drawcoastlines(linewidth=0.5, color='black')
    m.drawcountries(linewidth=0.5, color='black')

    x, y = m(lon, lat)

    m.contourf(x,y,vadv, cmap='jet',levels=16, \
               vmin=-2.e-8, vmax=2.e-8, extend='both')
    # Add a colorbar and title                                                                                      
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])         
    plt.colorbar(label='vorticity')
    plt.title('vorticity advection ' + time_str)
    plt.savefig("vort_adv"+time_str+".png")
    if showplot:
        plt.show()
    plt.close()

def plot_all_fields(dataset, nsteps, showplot):
    # plot the fields for all time periods
# Create a new map projection
    steps = dataset['step'].values
    for time_index in range(0, nsteps):

        steps = dataset['step'].values
        time_step = steps[time_index]
        print("time index is ", time_index)
        start_time = dataset['time'].values
        time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_index))
        time_str = np.datetime_as_string(time_stamp, unit='s')
        plot_vort_advection(data, time_index, time_str, showplot)
        plot_winds_v2(data,time_index, time_str, showplot) 
        plot_speed_v2(data,time_index, time_str, showplot)
        # plot the absolute  voriticity
        plot_vort_v2(data,time_index, time_str, showplot)
        # plot relative vorticity
        plot_rel_vort_v2(data,time_index, time_str, showplot)

# done making and saving plots.
#
# interpolate vorticity into each node of each trajectory
#
#def interp_vorticity(dataset, trajectory_list):
#
#
#    for trajectory in trajectory_list:
#        for p in trajectory.points:
#            p.vort = dataset['abs_vorticity'].interp(lat=p.lat, lon=p.lon,
#                                                 time=p.timestamp)
#            
#    print("Vorticities were added to trajectoriea")
#    return

def plot_vort_timeline(trajectories):
    """
        plot the vorticity as function of time for  each trajectory
    """

    max_points = 11 # limit on x axis - should match maximum number of points
    for trajectory in trajectories:
        if trajectory.length < 4:
            continue  # skip any empty or trunjcated edge  trajectories
        vort_list = []

        for index, point in enumerate(trajectory.points):
            if index == 0:
                start_lat = point.lat
                start_lon = point.lon
                start_time = point.timestamp
            stop_lat = point.lat
            stop_lon = point.lon
            stop_time = point.timestamp
            vort_list.append(point.vort)

        plt.plot(vort_list,marker="o")
        plt.xlabel('Time')
        plt.ylabel('Absolute Vorticity')
        plt.ylim(0, 2.e-4)
        plt.xlim(0,max_points)
        #
        # show start and stop lat and lon in title
        #
        title = "start:  {:.2f},{:.2f}, end: {:.2f},{:.2f}".\
            format(start_lat, start_lon, stop_lat, stop_lon)
        start_str = np.datetime_as_string(start_time, unit='s')
        stop_str =  np.datetime_as_string(stop_time, unit='s')
        title = title + "\n" + start_str + " to " + stop_str
        plt.title(title)

        # Display the plot
        plt.show()
        

    return # nothing 

# TODO: Analytics on changes in vorticity
# Try an RMS percentage change in the vorticity over all complete trajectorie
# Compare the vorticity at the start of a trajectory with its value at the
# end.
#
#    vort0 = []
#    vort1 = []
#    for trajectory in trajectory_list:
#        v0 = trajectory.points[0].vort
#        v1 = trajectory.points[-1].vort
#        if( math.isnan(v0)):
#            continue
#        if( math.isnan(v1)):
#            continue
#        if( v0 < 1.e-10):
#            continue
#        if( v1 < 1.e-10):
#            continue
#        vort0.append(v0)
#        vort1.append(v1)
#
#    
#    vort0 = np.array(vort0)
#    vort1 = np.array(vort1)
#    diff = (vort1 - vort0)**2
#    mean_diff = np.sum(diff) / len(vort1)
#    rms = np.sqrt(mean_diff)
#    print("RMS change in trajectory vorticity is ", rms)
#    rms0 = np.sqrt(np.sum(vort0* vort0)/len(vort0)        )
#    print("RMS of initial vorticity is ", rms0)
#    print("Relative change is ", rms/rms0*100, " percent")
def plot_traj_timeline(trajectory):
    fig, axs = plt.subplots(2, 1)  # 2 rows, 1 column

    # Plot on the first subplot
    #    axs[0].plot(x, y1)
    axs[0].set_title('Plot 1')

    # Plot on the second subplot
#    axs[1].plot(x, y2)
    axs[1].set_title('Plot 2')
    
    # Adjust spacing between subplots
    plt.tight_layout()

    # Display the subplots
    plt.show()

#
# interpolate the height data from the dataset at 1 deg lat
# to new one at 0.5 deg
def interp_z(ds, time, step, lats, lons):

    h = np.empty((len(lats), len(lons)))
    d = ds['z500'].sel(time=time, step=step)
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):

            z = d.interp(lat=lat, lon=lon, method='slinear')
    
            h[i][j] = z
                
                                        
    return h

# main program:
# Read winds and vorticity from NetCDF file if it exists. Else create it.
#
# repeat for each step up to 48 hours out.
maxsteps = 48  # maximum number of time steps we want to analyze

dataset_file = 'vortdata_steps.nc'
if os.path.exists(dataset_file):
    data = xr.open_dataset(dataset_file)
    print("Data read in from ", dataset_file)
else:
    # create the data and store it
    fc_ds_baro = xr.open_dataset("forecasts_baro-spharm_2007-04-15_2007-04-18.nc")
    fc_ds_baro = fc_ds_baro.roll({"lon": 180})
    steps =  fc_ds_baro.variables['step']
    nsteps = len(steps)

    #
    # create a new dataset with the winds, vorticity in it.
    timestamps = fc_ds_baro["time"].values
    start_time = timestamps[0]
    print("Start time: ", start_time)
    
    data = xr.Dataset(
        {
        "z500": (["step", "lat", "lon"], \
                  np.zeros((len(steps),  len(lat_lin), \
                            len(lon_lin)))), \

        "wind_u": (["step", "lat", "lon"], \
                  np.zeros(( len(steps), len(lat_lin), \
                            len(lon_lin)))), \
        "wind_v": (["step", "lat", "lon"], \
                  np.zeros(( len(steps),  len(lat_lin), \
                            len(lon_lin)))), \
        "rel_vorticity": (["step","lat", "lon"], \
                  np.zeros((len(steps), len(lat_lin), \
                                len(lon_lin)))),
        "abs_vorticity": (["step","lat", "lon"], \
                      np.zeros((len(steps),   len(lat_lin), \
                                len(lon_lin)))),

        "speed": (["step","lat", "lon"], \
                      np.zeros((len(steps),  len(lat_lin), \
                                len(lon_lin)))), \
        },\
        coords={"time":start_time, "step": ("step", steps), \
                "lat": ("lat",lat_lin), "lon": ("lon", lon_lin)}, \
    )

    # repeat for each time step up to 2 days
    # read in the data, compute winds and vorticity, add to dataset
# start with time step 1 - don't use 0, it is the initial data that is bad
    for step in range(0,maxsteps):
        time_step = steps[step].values
        print("time step is ", time_step)

        dt_str, time_stamp = get_timestamp(start_time, step)
        print("time stamp ", dt_str)
        # overwrite the height data with test data.

        # query the height data from the dataset
        geopot = baro_fcst(start_time, steps[step]) 
# interp_z is not working. dont use it.
#        geopot = interp_z(fc_ds_baro, start_time, steps[step],\
#                          lat_lin, lon_lin)
# these test patterns do work.
#        geopot = south_wind_v2(geopot, lat_lin, lon_lin)
#        geopot = north_wind_v2(geopot, lat_lin, lon_lin)
        # geopot = ridge_and_trough()
# verion 3 uses mu with gradient calls
#        winds_u, winds_v, zeta3, speed3 = prescribe_winds3()
        #        version 2 uses my centereed differences≈
        winds_u, winds_v, zeta3, speed3 = prescribe_winds2()
        #
        zeta = zeta3
        speed = speed3
        data['z500'][step] = geopot 
        data['wind_u'][step] = winds_u
        data['wind_v'][step] = winds_v
        data['rel_vorticity'][step] = zeta
        data['abs_vorticity'][step] = zeta + f
        data['speed'][step] = speed    

    data.to_netcdf(dataset_file)
    print("Data written to ", dataset_file)

m = Basemap(projection='cyl', llcrnrlat=min_lat, \
                urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon)

#plot_all_fields(data, maxsteps, False)

# Plotting the trajectories
#
print("All times completed")

trajeotory_file = "trajectories.nc"

start_time = data['time']

#

dt_hours = 1 # time step in hours
nsteps = 48 # number of step to integrate over

deltat = dt_hours * 3600 # time step in seconds
print("Computing trajectories")
trajectories = compute_trajectories(data, start_time, deltat, nsteps)
# plot trajectoriea at given time periods

print("Trajectory times:")
#
# trajectory end time is start time plus nsteps hours
start_time = trajectories[0].start_time
start_time_str = np.datetime_as_string(start_time, unit='s')
stop_time = np.datetime64(start_time + pd.Timedelta(hours=nsteps))
stop_time_str =  np.datetime_as_string(stop_time, unit='s')

print(start_time_str, stop_time_str)
#for i, t in enumerate(trajectories):
#    print(i, t.length, t.start_time, t.stop_time)
#    for p in t.points:
#        print("\t lat:", float(p.lat), " lon:", float(p.lon),
#              " dx:", float(p.dx), " dy: ", float(p.dy))
print("Plotting trajectories")
plot_trajectories(trajectories, start_time, stop_time)

#
plot_trajectories = True
if plot_trajectories:
    for i, t in enumerate(trajectories):
        filename = "tstep_" + str(dt_hours) + "hour" + str(i)
        plot_one_trajectory(t,filename, start_time, stop_time)
        print("saved ", filename)


