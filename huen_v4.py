import numpy as np
import pandas as pd
import copy
from traject_constants import *
from trajectorypoint import *
from trajectory_v2 import *
from velocity import Velocity
method='slinear' # interpolation method to use.
# chose from linear, slinear, quadratic, cubic
#
# get_timestamp - convert a starting time, and time_step in hours
def get_timestamp(start_time, time_step):
    time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_step))
    dt_str = np.datetime_as_string(time_stamp, unit='s')
    return dt_str, time_stamp


# get wind and vorticity interpolating in space linearly
def getwind_v4(ds, lat, lon, step):

    wind_u = ds['wind_u'][step].interp(lat=lat, lon=lon, method=method,
                                       assume_sorted=True)
    wind_v = ds['wind_v'][step].interp(lat=lat, lon=lon, method=method,
                                       assume_sorted=True)
    return copy.copy(wind_u), copy.copy(wind_v)
def getvort_v4(ds, lat, lon, step):
    vort   = ds['abs_vorticity'][step].interp(lat=lat, lon=lon, \
                                              method=method,
                                              assume_sorted=True)
    return copy.copy(vort)

#
# input the Xarray dataset, that can interpolated
def huen_v4(dataset, lat0, lon0, nsteps, deltat):
    """
    input;
         dataset
         lat0 - starting latitude for trajectory
         lon0 - starting longitude for trajectory
         nsteps - number of time periods to integrate over
         deltat - time step in seconds
    output:
         trajectory - list of trajectory nodes
    """
    #
    # in thie version, pass a iist of timestamps in netcdf format
    #
    # deltat is the time step in seconds between trajectory nodes.


    # lat0 and lon0 are the starting  point of the trajectory
    lat0:float; lat1:float; lon0: float; lon1: float
    dx:float; dy:float
    start_time_stamp = dataset['time'].values
    stop_time_str, stop_time_stamp = get_timestamp(start_time_stamp, nsteps)

    trajectory = Trajectory_v2(lat0, lon0, deltat, \
                               start_time_stamp, stop_time_stamp)
    lat1 = lat0
    lon1 = lon0  
# iterate through the timestamps, noting the last one terminates the
# the last node of the trajectory
    for timeindex in range(0,nsteps):
        next_time = timeindex + 1

#

# query the dstaset, interpolating in space
#
        wind_u, wind_v = getwind_v4(dataset, lat1, lon1, timeindex)

        # if no wind found at this location/time, break out of loop
        if np.isnan(wind_u) or np.isnan(wind_v): 
            break;

        denom = EARTH_RADIUS * np.cos(np.deg2rad(lat1))
        # lon_bar is the intermediate guess using euler method
        lon_bar  = lon1 + np.degrees(wind_u * deltat / denom)

        # lat_bar is intermediate guess using euler method
        lat_bar = lat1 + np.degrees(wind_v * deltat / EARTH_RADIUS )
        # get updated wind vector at updated longitude and latitude 

        wind_u_bar, wind_v_bar = \
            getwind_v4(dataset, lat_bar, lon_bar, next_time)
#        wind_u_bar, wind_v_bar = \
#            getwind(dataset, lat_bar, lon_bar, next_time)

        
        # again, if no wind found, break out of loop
        if np.isnan(wind_u_bar) or np.isnan(wind_v_bar) : 
            break;

        # update the longitude with huen's method.
        # use the euler updated latitude to convert distance to degrees               
        lon2 = lon1 + \
            0.5 * deltat *np.degrees((wind_u + wind_u_bar)/ \
                 (EARTH_RADIUS * np.cos(np.deg2rad(lat_bar))))

        lat2 = lat1 + \
            0.5 * deltat *np.degrees((wind_v + wind_v_bar) / EARTH_RADIUS)

        # calculate change in distance. Used for plotting trajectories
        # set to zero if in noise level
        dx = lon2 - lon1
        dy = lat2 - lat1
        if np.abs(dx) < 1.e-6:
            dx = 0
        if np.abs(dy) < 1.e-6:
            dy = 0
        # append the updated lat and lon to the trajectory object
        dt_str, time_stamp = get_timestamp(start_time_stamp, timeindex)

        point = TrajectoryPoint(lat1, lon1, dx, dy, timeindex)
        # interpolate to determine vorticity at this point.
        point.vort = getvort_v4(dataset, lat1, lon1, timeindex)
        # add the point to the trajectory
        trajectory.points.append(point)
        trajectory.length += 1  # keep track of length of trajectory
        #update the last position of the trajectory
        trajectory.last_lat = lat2
        trajectory.last_lon = lon2
        # advance to the next point
        lat1 = lat2
        lon1 = lon2 
#
# add the final point to the trajectory
#
    point = TrajectoryPoint(trajectory.last_lat, trajectory.last_lon,\
                            0., 0., timeindex)

    point.vort = getvort_v4(dataset,point.lat, point.lon, timeindex)

    trajectory.points.append(point)
    trajectory.length += 1
#
# update the trajectory stopping time using last point in trajectory
#
    t_str, trajectory.stop_time = get_timestamp(trajectory.start_time, timeindex)
#
## return the trajectory object
##                                                                                     
    return trajectory
    

        


        

                                   
