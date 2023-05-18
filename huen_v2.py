import numpy as np
from traject_constants import *
from trajectorypoint import *
from trajectory_v2 import *
from velocity import Velocity
        
    
def euler_v2(wind_model, lat0, lon0, deltat, nsteps):
    # same as huen_2 but with only the euler step included, for comparison
    trajectory = Trajectory_v2()
    # lat1 and lon1 are tracking the start point of the trajectory
    lat1 = lat0
    lon1 = lon0  

#
# repeat for each requested time step
##
    for time in range(0, nsteps):
        wind = wind_model(lat1, lon1, time)

        if wind is np.nan: # no wind found at this location/time
            break;

        denom = EARTH_RADIUS * np.cos(np.deg2rad(lat1))
        # lon_bar is the intermediate guess using euler method
        lon2  = lon1 + np.degrees(wind.u * deltat / denom)

        # lat_bar is intermediate guess using euler method
        lat2 = lat1 + np.degrees(wind.v * deltat / EARTH_RADIUS )
            
        dx = lon2 - lon1
        dy = lat2 - lat1  # used for plotting trajectories
        # append the updated lat and lon to the trajectory object
        point = TrajectoryPoint(lat1, lon1, dx, dy, time)
        trajectory.points.append(point)

        # advance to the next point
        lat1 = lat2
        lon1 = lon2 
#
## return the trajectory object
##                                                                                     
    return trajectory
#
# input the Xarray dataset, that can interpolated
def huen_v3(dataset, lat0, lon0, timestamps):
    """
    input;
         dataset
         lat0 - starting latitude for trajectory
         lon0 - starting longitude for trajectory
         timestamps - list of times for nodes in the trajectory
    output:
         trajectory - list of trajectory nodes
    """
    #
    # in thie version, pass a iist of timestamps in netcdf format
    #
    trajectory = Trajectory_v2()
    # lat1 and lon1 are tracking the start point of the trajectory
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
        deltat = int((next_time - timestamp)/ np.timedelta64(1, 's'))
#
# query the dstaset, interpolating in time and space
#
        wind_u = uwind.interp(time=timestamp, lat=lat1, \
                              lon=lon1, method="linear").values
        wind_v = vwind.interp(time=timestamp, lat=lat1, \
                              lon=lon1, method="linear").values

        # if no wind found at this location/time, break out of loop
        if np.isnan(wind_u) or np.isnan(wind_v): 
            break;

        denom = EARTH_RADIUS * np.cos(np.deg2rad(lat1))
        # lon_bar is the intermediate guess using euler method
        lon_bar  = lon1 + np.degrees(wind_u * deltat / denom)

        # lat_bar is intermediate guess using euler method
        lat_bar = lat1 + np.degrees(wind_v * deltat / EARTH_RADIUS )
        # get updated wind vector at updated longitude and latitude 

        wind_u_bar = uwind.interp(time=next_time, lat=lat_bar, \
                              lon=lon_bar, method="linear").values
        wind_v_bar = vwind.interp(time=next_time, lat=lat_bar, \
                              lon=lon_bar, method="linear").values

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
        trajectory.points.append(point)

        # advance to the next point
        lat1 = lat2
        lon1 = lon2 
#
## return the trajectory object
##                                                                                     
    return trajectory
    
def huen_v2(wind_model, lat0, lon0, deltat, nsteps):
    """
    huen - v2 - compute trajectory integrating by huens method                             
    input:
	wind_model - method that returns the velocity (m/s)
        at specified           
        input latitude, longitude and time step.
	This function is variable, to accomodate 
	test functions as well as for production wind data.                               
    input lat0 - starting latitude of parcel                                          
    input lon0 - starting longitude of parcel                                         
    input deltat - time step in seconds.                                              
    input nsteps - number of time steps to integrate over                                 
    Returns:                                                                          
    Trajectory_v2 object            
    """

    trajectory = Trajectory_v2()
    # lat1 and lon1 are tracking the start point of the trajectory
    lat1 = lat0
    lon1 = lon0  

#
# repeat for each requested time step
# adding in a starting point
#
    for time in range(0, nsteps +1):
        wind = wind_model(lat1, lon1, time)

        if wind is np.nan: # no wind found at this location/time
            break;

        denom = EARTH_RADIUS * np.cos(np.deg2rad(lat1))
        # lon_bar is the intermediate guess using euler method
        lon_bar  = lon1 + np.degrees(wind.u * deltat / denom)

        # lat_bar is intermediate guess using euler method
        lat_bar = lat1 + np.degrees(wind.v * deltat / EARTH_RADIUS )
        # get updated wind vector at updated longitude and latitude 
        # and at the next time step.                                
        wind_bar = wind_model(lat_bar, lon_bar, time + 1)

        if wind_bar is None: #  lat & lon are out of bounds          
            break;

        # update the longitude with huen's method.
        # use the euler updated latitude to convert distance to degrees               
        lon2 = lon1 + \
            0.5 * deltat *np.degrees((wind.u + wind_bar.u)/ \
                 (EARTH_RADIUS * np.cos(np.deg2rad(lat_bar))))

        lat2 = lat1 + \
            0.5 * deltat *np.degrees((wind.v + wind_bar.v) / EARTH_RADIUS)

            
        dx = lon2 - lon1
        dy = lat2 - lat1  # used for plotting trajectories
        # append the updated lat and lon to the trajectory object
        point = TrajectoryPoint(lat1, lon1, dx, dy, time)
        trajectory.points.append(point)

        # advance to the next point
        lat1 = lat2
        lon1 = lon2 
#
## return the trajectory object
##                                                                                     
    return trajectory

        


        

                                   
