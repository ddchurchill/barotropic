import numpy as np
import matplotlib.pyplot as plt
#
earth_radius = 6371e3 # in meters
    
# define constant geostrophic wind for testing.
# this method gets passed as an argument to the integration routine
# 
def geowind(lat, lon, time_step):
    return 10,1

#
# step_forward - compute trajectory by simple forward step
# input lat0 - starting latitude of parcel
# input lon0 - starting longitude of parcel
# input deltat - time step in seconds. 
# input nt - number of time steps to integrate over
# input velocity - the function that returns the u and v components of the wind
#
# output: lon_traj, lat_traj - linear arrays of length nt steps showing position of parcel
#
def step_forward(velocity, lat0, lon0, deltat, nt):
#
# initialize the trajectory
#
    lat_traj = np.zeros(nt, dtype=float)
    lon_traj = np.zeros(nt, dtype=float)
    lon_traj[0] = lon0
    lat_traj[0] = lat0
    for t in range(0, nt -1): # repeat for each time step

#        ug = getu(lat_traj[t], lon_traj[t], t)
#        vg = getv(lat_traj[t], lon_traj[t], t)
        ug, vg = velocity(lat_traj[t], lon_traj[t], t)
        dlambda = ug * deltat / (earth_radius * np.cos(lat_traj[t] * np.pi/180))* 180./np.pi
#        print("ug, lat, deltat, dlambda: ", ug, lat_traj[t], deltat, dlambda)
        dphi = vg * deltat /(earth_radius) * 180/np.pi
        lat_traj[t+1] = lat_traj[t] + dphi
        lon_traj[t+1] = lon_traj[t] + dlambda
        
    
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
# initialize the trajectory
#
    lat_traj = np.zeros(nt, dtype=float)
    lon_traj = np.zeros(nt, dtype=float)
    lon_traj[0] = lon0 # initialize the trajectory longitude
    lat_traj[0] = lat0  # initalize the trajectory latitude
    traj.lon =[lon0]
    traj.lat = [lat0]
#
# repeat for each time step
# by exit loop if the velocity was not found -- as when the lat and lon
# are out of bounds 
    for t in range(1, nt): 

        # get the wind components at the current lat, lon, and time step
#        u, v = velocity(lat_traj[t], lon_traj[t], t)
        u, v = velocity(traj.lat[-1], traj.lon[-1], t)

        # update the longitude
        #deltax = earth_radius * np.cos(np.deg2rad(lat_traj[t]))
        deltax = earth_radius * np.cos(np.deg2rad(traj.lat[-1]))
        # euler forward step first for longitude
#        lon_euler = lon_traj[t] + np.degrees(u * deltat / deltax)
        euler.lon = traj.lon[-1] + np.degrees(u * deltat / deltax)
        # get euler's update to latitude
#        lat_euler = lat_traj[t] + np.degrees(v * deltat / earth_radius )
        euler.lat = traj.lat[-1] + np.degrees(v * deltat / earth_radius )

        # get updated wind vector at updated longitude and latitude
        # and at the next time step.
        
        u_next, v_next = velocity(euler.lat, euler.lon, t+1)

        # update the longitude with huen's method.
        # use the euler updated latitude to convert distance to degrees
        deltax = earth_radius * np.cos(np.deg2rad(euler.lat))
        #        lon_traj[t+1] = lon_traj[t] \
        euler.lon = traj.lon[-1] + 0.5 * deltat *np.degrees((u + u_next)/deltax)
        
        # update the latitude with euler's method
        # lat_traj[t+1] = lat_traj[t] + \

        euler.lat = traj.lat[-1] + \
            0.5 * deltat *np.degrees((v + v_next) / earth_radius)

        # append the euler updated lat and lon to the arrays
        traj.lon.append(euler.lon)
        traj.lat.append(euler.lat)

    return traj.lat, traj.lon


def main():
    nt = 10
    deltat = 3600  # 1 hour integration period
    lat0 = 10.1 # degrees
    lon0 = 20.1 # degrees starting point
    lat_traj, lon_traj = step_forward(geowind, lat0, lon0, deltat, nt)
#    lat_traj, lon_traj = huen(lat0, lon0, deltat, nt)
    np.set_printoptions(precision=6)
    print("latitudes: ",lat_traj)
    print("longitudes: ", lon_traj)
    plt.figure()
    plt.plot(lon_traj, lat_traj)
    plt.show()
    
if __name__ == "__main__":
    main()          
