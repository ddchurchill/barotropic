import numpy as np
import matplotlib.pyplot as plt
from velocity import Velocity as vel
from trajectory import Trajectory 
#
earth_radius = 6371e3 # in meters
    
# define constant geostrophic wind for testing.
# this method gets passed as an argument to the integration routine
# 
def swwind(lat, lon, time_step):
    return vel(10,10)
def westwind(lat,lon, time_step):
    return vel(20, 0)
def eastwind(lat, lon, time_step):
    return vel(-20, 0)
def newind(lat, lon, time_step):
    return vel(-10, -10)
def northwind(lat, lon, time_step):
    return vel( 0, -10)

def cyclone(lat, lon, time_step):
    if(time_step == 0):
        return vel(20, 0)
    if(time_step == 1):
        return vel(20, 0)
    if(time_step == 2):
        return vel(10, 10)
    if(time_step == 3):
        return vel(10, 10)
    if(time_step == 4):
        return vel(0, 20)
    if(time_step == 5):
        return vel(0, 20)
    if(time_step == 6):
        return vel(-20, 0)
    if(time_step == 7):
        return vel(-20, 0)
    if(time_step == 8):
        return vel(-10, -10)
    if(time_step == 9):
        return vel(-10, -10)
    if(time_step == 10):
        return vel(0, -20)
    return vel(0,-20)

def anticyclone(lat, lon, time_step):
    if(time_step == 0):
        v = vel(-20,0)
        return v
    if(time_step == 1):
        return vel(-20, 0)
    if(time_step == 2):
        return vel(-10, -10)
    if(time_step == 3):
        return vel(-10, -10)
    if(time_step == 4):
        return vel(0, -20)
    if(time_step == 5):
        return vel(0, -20)
    if(time_step == 6):
        return vel(20, 0)
    if(time_step == 7):
        return vel(20, 0)
    if(time_step == 8):
        return vel(10, 10)
    if(time_step == 9):
        return vel(10, 10)
    if(time_step == 10):
        return vel(0, 20)
    return vel(0,20)

        

    
#
# nowind returns null values for winds, as when lat and lon are out of bounds
#
def  nowind(lat, lon, time_step):
    return None



def euler( velocity, lat0, lon0, deltat, nt):
    """
    euler - compute trajectory by simple forward step
    Parameters:
    input lat0 - starting latitude of parcel
    input lon0 - starting longitude of parcel
     input deltat - time step in seconds. 
    input nt - number of time steps to integrate over
    input velocity - the function that returns the u and v components of the wind
    at the point specified by the argument.

    Returns:
    Trajectory object - linear arrays of length nt steps showing position of parcel
    """
#
#
# initialize the trajectory
#
    traj = Trajectory(lat0, lon0)
    for t in range(0, nt): # repeat for each time step

        wind = velocity(traj.last_lat(), traj.last_lon(), t)
        if wind is None:
            break;
        deltax = earth_radius * np.cos(np.deg2rad(traj.last_lat()))


        dphi = np.degrees(wind.v * deltat /earth_radius) 
        euler_lat = traj.last_lat() + dphi
        traj.lat.append( euler_lat)
        
        dlambda = wind.u * deltat / deltax
        euler_lon = traj.last_lon() + np.degrees( dlambda)
        traj.lon.append(euler_lon)
        
    
    return traj

#
def huen( wind_func, lat0, lon0, deltat, nt):
    """ 
    huen - compute trajectory integrating by huens method
    input wind_func - method that returns the velocity at specified
     input latitude, longitude and time step. This function is variable, to accomodate
    test functions as well as for production wind data. 
    input lat0 - starting latitude of parcel
    input lon0 - starting longitude of parcel
    input deltat - time step in seconds. 
    input nt - number of time steps to integrate over

    Returns:
    Trajectory object
"""
#
# initialize the trajectory with the input lat and lon
#

    traj = Trajectory(lat0, lon0)
#
# repeat for each time step
# but exit loop if the velocity was not found -- as when the lat and lon
# are out of bounds 
    for t in range(0, nt +1): 

        # get the wind components at the current lat, lon, and time step
#        print("huen lat, lon, time:", traj.last_lat(), traj.last_lon(), t)
        vector = wind_func(traj.last_lat(), traj.last_lon(), t)
        if vector is np.nan:
            break
        u = vector.u
        v = vector.v

        if u is None:
            break;

        # update the longitude
        deltax = earth_radius * np.cos(np.deg2rad(traj.last_lat()))
        # euler forward step first for longitude
        euler_lon = traj.last_lon() + np.degrees(u * deltat / deltax)
        # get euler's update to latitude
        euler_lat = traj.last_lat() + np.degrees(v * deltat / earth_radius )

        # get updated wind vector at updated longitude and latitude
        # and at the next time step.
        
        next_wind = wind_func(euler_lat, euler_lon, t+1)

        if next_wind is None: # no wind found -- lat & lon are out of bounds
            break;

        # update the longitude with huen's method.
        # use the euler updated latitude to convert distance to degrees
        deltax = earth_radius * np.cos(np.deg2rad(euler_lat))

        euler_lon = traj.last_lon() + 0.5 * deltat *np.degrees((u + next_wind.u)/deltax)
        
        # update the latitude with euler's method

        euler_lat = traj.last_lat() + \
            0.5 * deltat *np.degrees((v + next_wind.v) / earth_radius)

        # append the euler-updated lat and lon to the trajectory object
        traj.lon.append(euler_lon)
        traj.lat.append(euler_lat)
        traj.length = traj.length + 1
#
# return the trajectory object
#
    return traj


def main():
    nt = 12
    deltat = 3600  # 1 hour integration period
    lat0 = 10. # degrees
    lon0 = 20. # degrees starting point

    traj = huen(anticyclone, lat0, lon0, deltat, nt)

#    traj = euler(anticyclone, lat0, lon0, deltat, nt)

#    lat_traj, lon_traj = euler(cyclone, lat0, lon0, deltat, nt)
    np.set_printoptions(precision=6)

# plot the trajectory
    for i in range(len(traj.lat)-1):
        dx = traj.lon[i+1] - traj.lon[i]
        dy = traj.lat[i+1] - traj.lat[i]
        plt.arrow(traj.lon[i], traj.lat[i], dx, dy, \
                  length_includes_head=True, head_length=0.1, head_width=0.05)

    
    plt.show()
#
# now test when no winds are found
#
#    lat_traj, lon_traj = huen(nowind, lat0, lon0, deltat, nt)
#    print("no wind lat: ",lat_traj)
#    print("no wind lon: ", lon_traj)
    
if __name__ == "__main__":
    main()          
