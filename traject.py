import numpy as np
import matplotlib.pyplot as plt
#
earth_radius = 6371e3 # in meters
    
# define constant geostrophic wind for testing.
# this method gets passed as an argument to the integration routine
# 
def swwind(lat, lon, time_step):
    return 10,10
def westwind(lat,lon, time_step):
    return 20, 0
def eastwind(lat, lon, time_step):
    return -20, 0
def newind(lat, lon, time_step):
    return -10, -10
def northwind(lat, lon, time_step):
    return 0, -10

def cyclone(lat, lon, time_step):
    if(time_step == 0):
        return 20, 0
    if(time_step == 1):
        return 20, 0
    if(time_step == 2):
        return 10, 10
    if(time_step == 3):
        return 10, 10
    if(time_step == 4):
        return 0, 20
    if(time_step == 5):
        return 0, 20
    if(time_step == 6):
        return -20, 0
    if(time_step == 7):
        return -20, 0
    if(time_step == 8):
        return -10, -10
    if(time_step == 9):
        return -10, -10
    if(time_step == 10):
        return 0, -20
    return 0,-20

def anticyclone(lat, lon, time_step):
    if(time_step == 0):
        return -20, 0
    if(time_step == 1):
        return -20, 0
    if(time_step == 2):
        return -10, -10
    if(time_step == 3):
        return -10, -10
    if(time_step == 4):
        return 0, -20
    if(time_step == 5):
        return 0, -20
    if(time_step == 6):
        return 20, 0
    if(time_step == 7):
        return 20, 0
    if(time_step == 8):
        return 10, 10
    if(time_step == 9):
        return 10, 10
    if(time_step == 10):
        return 0, 20
    return 0,20

        

    
#
# nowind returns null values for winds, as when lat and lon are out of bounds
#
def nowind(lat, lon, time_step):
    return None, None
#
# euler - compute trajectory by simple forward step
# input lat0 - starting latitude of parcel
# input lon0 - starting longitude of parcel
# input deltat - time step in seconds. 
# input nt - number of time steps to integrate over
# input velocity - the function that returns the u and v components of the wind
#
# output: lon_traj, lat_traj - linear arrays of length nt steps showing position of parcel
#
def euler(velocity, lat0, lon0, deltat, nt):
#
# initialize the trajectory
#
#    lat_traj = np.zeros(nt, dtype=float)
#    lon_traj = np.zeros(nt, dtype=float)
    lon_traj = [lon0]
    lat_traj = [lat0]
    for t in range(0, nt): # repeat for each time step

        ug, vg = velocity(lat_traj[-1], lon_traj[-1], t)
        if ug is None:
            break;
        deltax = earth_radius * np.cos(np.deg2rad(lat_traj[-1]))
        dlambda = ug * deltat / deltax

        dphi = np.degrees(vg * deltat /earth_radius) 
        euler_lat = lat_traj[-1] + dphi
        lat_traj.append( euler_lat)
        euler_lon = lon_traj[-1] + np.degrees( dlambda)
        lon_traj.append(euler_lon)
        
    
    return lat_traj, lon_traj

# huen - compute trajectory integrating by huens method
# input velocity - method that returns the velocity at specified
#     input latitude, longitude and time step. This function is variable, to accomodate
# test functions as well as for production wind data. 
# input lat0 - starting latitude of parcel
# input lon0 - starting longitude of parcel
# input deltat - time step in seconds. 
# input nt - number of time steps to integrate over
#
# output: lon_traj, lat_traj - linear arrays of length nt steps showing position of parcel
#
def huen(velocity, lat0, lon0, deltat, nt):

#
# initialize the trajectory with the input lat and lon
#

    traj_lon = [lon0]
    traj_lat = [lat0]
#
# repeat for each time step
# but exit loop if the velocity was not found -- as when the lat and lon
# are out of bounds 
    for t in range(0, nt): 

        # get the wind components at the current lat, lon, and time step
        u, v = velocity(traj_lat[-1], traj_lon[-1], t)

        if u is None:
            break;

        # update the longitude
        deltax = earth_radius * np.cos(np.deg2rad(traj_lat[-1]))
        # euler forward step first for longitude
        euler_lon = traj_lon[-1] + np.degrees(u * deltat / deltax)
        # get euler's update to latitude
        euler_lat = traj_lat[-1] + np.degrees(v * deltat / earth_radius )

        # get updated wind vector at updated longitude and latitude
        # and at the next time step.
        
        u_next, v_next = velocity(euler_lat, euler_lon, t+1)

        if u_next is None: # no velocity found -- lat & lon are out of bounds
            break;

        # update the longitude with huen's method.
        # use the euler updated latitude to convert distance to degrees
        deltax = earth_radius * np.cos(np.deg2rad(euler_lat))

        euler_lon = traj_lon[-1] + 0.5 * deltat *np.degrees((u + u_next)/deltax)
        
        # update the latitude with euler's method

        euler_lat = traj_lat[-1] + \
            0.5 * deltat *np.degrees((v + v_next) / earth_radius)

        # append the euler updated lat and lon to the arrays
        traj_lon.append(euler_lon)
        traj_lat.append(euler_lat)
#
# return the arrays
#
    return traj_lat, traj_lon


def main():
    nt = 12
    deltat = 3600  # 1 hour integration period
    lat0 = 10. # degrees
    lon0 = 20. # degrees starting point
#    lat_traj, lon_traj = huen(newind, lat0, lon0, deltat, nt)
#    lat_traj, lon_traj = huen(cyclone, lat0, lon0, deltat, nt)
    lat_traj, lon_traj = huen(anticyclone, lat0, lon0, deltat, nt)
#    lat_traj, lon_traj = euler(anticyclone, lat0, lon0, deltat, nt)
#    lat_traj, lon_traj = euler(cyclone, lat0, lon0, deltat, nt)
    np.set_printoptions(precision=6)

# plot the trajectory
    for i in range(len(lat_traj)-1):
        dx = lon_traj[i+1] - lon_traj[i]
        dy = lat_traj[i+1] - lat_traj[i]
        plt.arrow(lon_traj[i], lat_traj[i], lon_traj[i+1] -lon_traj[i], lat_traj[i+1]-lat_traj[i], \
                  length_includes_head=True, head_length=0.1, head_width=0.05)

#    plt.figure()
#    plt.plot(lon_traj, lat_traj)
    
    plt.show()
#
# now test when no winds are found
#
#    lat_traj, lon_traj = huen(nowind, lat0, lon0, deltat, nt)
#    print("no wind lat: ",lat_traj)
#    print("no wind lon: ", lon_traj)
    
if __name__ == "__main__":
    main()          
