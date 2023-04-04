import numpy as np
#
# get geostrophic wind returns the geostrophic wind at given latitude, longitude, and time_step
#
# returns:  ug, vg - components of the geostropic wind
#
def get_geostrophic_wind(lat, lon, time_step):
    ug = 100
    vg = 0
    # TODO: interpolate against the grided wind to get the precise
    # components
    return ug, vg


#
# compute trajectory by huens method
# input lat0 - starting latitude of parcel
# input lon0 - starting longitude of parcel
# input deltat - time step in seconds. 
# input nt - number of time steps to integrate over
#
# output: lon_traj, lat_traj - linear arrays of length nt steps showing position of parcel
#
def integrate(lat0, lon0, deltat, nt):
    earth_radius = 6371e3 # in meters
#
# initialize the trajectory
#
    lat_traj = np.zeros(nt, dtype=float)
    lon_traj = np.zeros(nt, dtype=float)
    lon_traj[0] = lon0
    lat_traj[0] = lat0
    for t in range(0, nt -1): # repeat for each time step

        ug, vg = get_geostrophic_wind(lat_traj[t], lon_traj[t], t)

        dlambda = ug * deltat / (earth_radius * np.cos(lat_traj[t] * np.pi/180))* 180./np.pi
#        dlambda = ug * deltat / 111e3 # 111 km per degree of longitude roughly
        print("ug, lat, deltat, dlambda: ", ug, lat_traj[t], deltat, dlambda)
        dphi = vg * deltat /(earth_radius) * 180/np.pi
        lat_traj[t+1] = lat_traj[t] + dphi
        lon_traj[t+1] = lon_traj[t] + dlambda
        
    
    return lat_traj, lon_traj

def main():
    nt = 10
    deltat = 3600  # 1 hour integration period
    lat0 = 10.1 # degrees
    lon0 = 20.1 # degrees starting point
    lat_traj, lon_traj = integrate(lat0, lon0, deltat, nt)
    np.set_printoptions(precision=6)
    print("latitudes: ",lat_traj)
    print("longitudes: ", lon_traj)
    
if __name__ == "__main__":
    main()          
