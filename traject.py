import numpy as np
import matplotlib.pyplot as plt
#
earth_radius = 6371e3 # in meters
    
# get geostrophic wind returns the geostrophic wind at given latitude, longitude, and time_step
#
# returns:  ug, vg - components of the geostropic wind
#
#def get_geostrophic_wind(lat, lon, time_step):
#    ug = 100
#    vg = 50
#    # TODO: interpolate against the grided wind to get the precise
#    # components
#    return ug, vg
def getu(lat, lon, time_step): # get u component of geostrophic wind
#    # TODO: interpolate against the grided wind to get the precise
    return 10
def getv(lat, lon, time_step): # get v component of geostrophic wind
    # TODO: interpolate against the grided wind to get the precise wind
    return 10

#
# step_forward - compute trajectory by simple forward step
# input lat0 - starting latitude of parcel
# input lon0 - starting longitude of parcel
# input deltat - time step in seconds. 
# input nt - number of time steps to integrate over
#
# output: lon_traj, lat_traj - linear arrays of length nt steps showing position of parcel
#
def step_forward(lat0, lon0, deltat, nt):
#
# initialize the trajectory
#
    lat_traj = np.zeros(nt, dtype=float)
    lon_traj = np.zeros(nt, dtype=float)
    lon_traj[0] = lon0
    lat_traj[0] = lat0
    for t in range(0, nt -1): # repeat for each time step

        ug = getu(lat_traj[t], lon_traj[t], t)
        vg = getv(lat_traj[t], lon_traj[t], t)
        dlambda = ug * deltat / (earth_radius * np.cos(lat_traj[t] * np.pi/180))* 180./np.pi
#        print("ug, lat, deltat, dlambda: ", ug, lat_traj[t], deltat, dlambda)
        dphi = vg * deltat /(earth_radius) * 180/np.pi
        lat_traj[t+1] = lat_traj[t] + dphi
        lon_traj[t+1] = lon_traj[t] + dlambda
        
    
    return lat_traj, lon_traj

# huen - compute trajectory integrating by huens method
# input lat0 - starting latitude of parcel
# input lon0 - starting longitude of parcel
# input deltat - time step in seconds. 
# input nt - number of time steps to integrate over
#
# output: lon_traj, lat_traj - linear arrays of length nt steps showing position of parcel
#
def huen(lat0, lon0, deltat, nt):

#
# initialize the trajectory
#
    lat_traj = np.zeros(nt, dtype=float)
    lon_traj = np.zeros(nt, dtype=float)
    lon_traj[0] = lon0
    lat_traj[0] = lat0
    for t in range(0, nt -1): # repeat for each time step

        ug = getu(lat_traj[t], lon_traj[t], t)
        deltax = earth_radius * np.cos(np.deg2rad(lat_traj[t]))
        lon_bar = lon_traj[t] + np.degrees(ug * deltat / deltax)
        lon_traj[t+1] = lon_traj[t] \
            + 0.5 * deltat *np.degrees((ug + getu(lat_traj[t], lon_bar, t+1))/deltax)                                                 
        vg = getv(lat_traj[t], lon_traj[t], t)
        lat_bar = lat_traj[t] + np.degrees(vg * deltat / earth_radius )
        lat_traj[t+1] = lat_traj[t] + \
            0.5 * deltat *np.degrees((vg + getv(lat_bar,lon_traj[t], t +1)) / earth_radius)
    return lat_traj, lon_traj


def main():
    nt = 10
    deltat = 3600  # 1 hour integration period
    lat0 = 10.1 # degrees
    lon0 = 20.1 # degrees starting point
#    lat_traj, lon_traj = step_forward(lat0, lon0, deltat, nt)
    lat_traj, lon_traj = huen(lat0, lon0, deltat, nt)
    np.set_printoptions(precision=6)
    print("latitudes: ",lat_traj)
    print("longitudes: ", lon_traj)
    plt.figure()
    plt.plot(lon_traj, lat_traj)
    plt.show()
    
if __name__ == "__main__":
    main()          
