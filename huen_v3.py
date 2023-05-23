import numpy as np
from traject_constants import *
from trajectorypoint import *
from trajectory_v2 import *
from velocity import Velocity

def getwind(ds, lat, lon, timestamp):
    wind_u = ds['wind_u'].interp(time=timestamp, lat=lat, \
                          lon=lon, method="linear").values
    wind_v = ds['wind_v'].interp(time=timestamp, lat=lat, \
                          lon=lon, method="linear").values
    return wind_u, wind_v

def getvort(ds, lat, lon, timestamp):
    vort = ds['abs_vorticity'].interp(lat=lat, \
                        lon=lon, time=timestamp)
    return vort


    
    
#
# input the Xarray dataset, that can interpolated
def huen_v3(dataset, lat0, lon0, timestamps, deltat):
    """
    input;
         dataset
         lat0 - starting latitude for trajectory
         lon0 - starting longitude for trajectory
         timestamps - list of times for nodes in the trajectory
         deltat - time step in seconds
    output:
         trajectory - list of trajectory nodes
    """
    #
    # in thie version, pass a iist of timestamps in netcdf format
    #
    # deltat is the time step in seconds between trajectory nodes.
#    deltat = int((timestamps[1] - timestamps[0])/ np.timedelta64(1, 's'))

    # lat0 and lon0 are the starting  point of the trajectory
    trajectory = Trajectory_v2(lat0, lon0, deltat, timestamps[0])

    lat1 = lat0
    lon1 = lon0  
# iterate through the timestamps, noting the last one terminates the
# the last node of the trajectory
    uwind = dataset['wind_u']
    vwind = dataset['wind_v']
    for timeindex, timestamp in enumerate(timestamps[:-1]):
        next_time = timestamps[timeindex +1]
#
# compute elapsed time in seconds between requested timestamps
#
#
# query the dstaset, interpolating in time and space
#
        wind_u, wind_v = getwind(dataset, lat1, lon1, timestamp)
        # if no wind found at this location/time, break out of loop
        if np.isnan(wind_u) or np.isnan(wind_v): 
            break;

        denom = EARTH_RADIUS * np.cos(np.deg2rad(lat1))
        # lon_bar is the intermediate guess using euler method
        lon_bar  = lon1 + np.degrees(wind_u * deltat / denom)

        # lat_bar is intermediate guess using euler method
        lat_bar = lat1 + np.degrees(wind_v * deltat / EARTH_RADIUS )
        # get updated wind vector at updated longitude and latitude 

#        wind_u_bar = uwind.interp(time=next_time, lat=lat_bar, \
#                              lon=lon_bar, method="linear").values
#        wind_v_bar = vwind.interp(time=next_time, lat=lat_bar, \
#                              lon=lon_bar, method="linear").values
        wind_u_bar, wind_v_bar = getwind(dataset, lat_bar, lon_bar, next_time)
        
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

            
        dx = lon2 - lon1
        dy = lat2 - lat1  # used for plotting trajectories
        # append the updated lat and lon to the trajectory object
        point = TrajectoryPoint(lat1, lon1, dx, dy, timestamp)
        # interpolate to determine vorticity at this point.
#        point.vort = dataset['abs_vorticity'].interp(lat=point.lat, \
#                                                     lon=point.lon, \
#                                                 time=point.timestamp)
        point.vort = getvort(dataset, lat1, lon1, timestamp)
        # add the point to the trajectory
        trajectory.points.append(point)
        trajectory.length += 1  # keep track of length of trajectory
        #update the last position of the trajectory
        trajectory.last_lat = lat2
        trajectory.last_lon = lon2
        trajectory.stop_time = next_time
        # advance to the next point
        lat1 = lat2
        lon1 = lon2 
#
# add the final point to the trajectory
#
    point = TrajectoryPoint(trajectory.last_lat, trajectory.last_lon,\
                            0., 0., trajectory.stop_time)
#    point.vort = dataset['abs_vorticity'].interp(lat=point.lat, \
#                                                     lon=point.lon, \
#                                                 time=point.timestamp)
    point.vort = getvort(dataset,point.lat, point.lon, point.timestamp)
    trajectory.points.append(point)
    trajectory.length += 1
#
## return the trajectory object
##                                                                                     
    return trajectory
    

        


        

                                   
