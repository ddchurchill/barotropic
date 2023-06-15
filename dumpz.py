import pprint
import xarray as xr
import numpy as np
import pandas as pd
from traject_constants import * 
def laplacian(data, x, y):
    nrows, ncols = data.shape
    laplace = np.zeros((nrows, ncols))
    delsquared = np.zeros((nrows,ncols))
    dx = np.zeros((nrows, ncols))
    ddx = np.zeros((nrows, ncols))
    d2x = np.zeros((nrows, ncols))
    d2y = np.zeros((nrows, ncols))    

    dy = np.zeros((nrows, ncols))
    ddy = np.zeros((nrows, ncols))    
    print("size: ", nrows, ncols)
    for i in 1, nrows -2:
        for j in 1, ncols -2:
            dx[i,j] = x[i,j+1] - x[i,j]
            dy[i,j] = y[i+1,j] - y[i,j]
            ddx[i,j] = (data[i,j+1] - data[i,j-1])/(2 * dx[i,j])
            ddy[i,j] = (data[i+1,j] - data[i-1,j])/(2 * dy[i,j])
            d2x[i,j] =  data[i,j+1] - 2*data[i,j] + data[i,j-1]
            d2y[i,j] = data[i+1,j] - 2*data[i,j] + data[i-1,j]
            delsquared[i,j] = d2x[i,j] + d2y[i,j]
            laplace[i,j] = d2x[i,j]/dx[i,j]**2 + d2y[i,j]/dy[i,j]**2            

    return laplace, dx, dy, ddx, ddy, d2x, d2y, delsquared

def get_timestamp(start_time, time_step):
    time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_step))
    dt_str = np.datetime_as_string(time_stamp, unit='s')
    return dt_str, time_stamp

def print_array(title, data):
    line = "_" * 60
    print()

    print("\t\t", title)
    print(line)
    print("\t\t\tLongitude")
    print("Lat ", end="")
    for i in range(0,nlon):
        print(f"\t{lon_lin[i]:10.2f}", end="")
    print("")
    print("-"*60)
    for i in range(nlat -1, -1, -1):
            print(lat_lin[i], end="")
            for j in range(0, nlon):
                d = data[i,j]
                print(f"\t{d:.4e}", end=" ")
            print("")

    print("-"*60)


#print("lat_lin: ",lat_lin)
if __name__ == "__main__":

    min_lon = -56
    max_lon = -54
    min_lat = 34
    max_lat = 36
    nlat = max_lat - min_lat +1
    
    nlon = max_lon - min_lon +1
    print(nlat, nlon)
    lat_lin = np.linspace(min_lat, max_lat, nlat)
    lon_lin = np.linspace(min_lon, max_lon, nlon)
    lon, lat = np.meshgrid(lon_lin, lat_lin)
    phi_rad = np.array(lat_lin * np.pi/180.)  # latitude in radians
    lambda_rad = np.array(lon_lin * np.pi/180.)  # longitude in radians
    lat_rad = np.radians(lat)  # convert degrees to radians
    lon_rad = np.radians(lon)
    f = 2 * OMEGA * np.sin(lat_rad)  # Coriolis parameter
    x = np.cos(np.deg2rad(lat))*EARTH_RADIUS*lon_rad
    y = EARTH_RADIUS * np.deg2rad(lat)
    #
    print_array('X', x)
    print_array('Y', y)#
    
    dataset = xr.open_dataset("forecasts_baro-spharm_2007-04-15_2007-04-18.nc")
#    dataset = dataset.roll({"lon": 180})
    
    steps =  dataset.variables['step']
#    print("steps", steps)
#    nsteps = len(steps)
    step = 7 # 7z on 15th 
    
#    timestamps = dataset["time"].values.copy()
#    forecast_init_time = timestamps[0] # 15Sep 2007
#    start_time_str, start_time = get_timestamp(forecast_init_time, step)
#    print("start time: ", start_time_str)
    #
    # get a window of height data, centered over a jet streak
    #
    #z500 = dataset.sel({
    #    "time": forecast_init_time,
    #    "step": steps[step],
    #    "lat": slice(max_lat, min_lat),
    #    "lon": slice(180+min_lon, 180+max_lon)
    #})['z500'].values.copy()
    
    #z500 = np.flip(z500,axis=0) # flip the data as in geowinds4.
    #
    # now get the derived data set
    #
    derived = xr.open_dataset('vortdata_steps.nc')
    
    h500 = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['z500'].values.copy()
    
    looked_rel_vorticity = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['rel_vorticity'].values.copy()
    
    #
    # h500 and z500 should be identical
    #
    #mse = np.mean((h500 - z500) ** 2)
    #rmse = np.sqrt(mse)
    #if rmse != 0:
    #    print("Error: z500 and h500 differ")
    #    sys.exit()
    #
    # print a tables of data
    #
    #print_array('z500', z500)
    print_array('h500', h500)
    
    
    zeta, dx, dy, ddx, ddy, d2x, d2y, delsquared  = laplacian(h500, x, y)
    
    print_array("Dx", dx)
    print_array('Dy', dy)
    print_array('ddx', ddx)
    print_array('ddy', ddy)
    print_array('d2x', d2x)
    print_array('d2y', d2y)
    print_array('Laplacian', zeta)
    print_array('Del squared', delsquared)
    print_array('coriolis', f)
    
    # compute vorticity, multiply by gravity, divide by f
    computed_rel_vorticity = zeta * GRAVITY / f
    print_array('Computed Vorticity', computed_rel_vorticity)
    f_35deg = 2. * OMEGA * np.sin(np.deg2rad(35.))
    print("Omega: ", OMEGA)
    print("f 35 deg: ", f_35deg)
    print_array('Looked up Rel vorticity', looked_rel_vorticity)
    
    u = -GRAVITY/f*ddy
    v = GRAVITY/f*ddx
    
    speed = np.sqrt(u*u + v*v)
    print_array('Computed U', u)
    print_array('Computed V', v)
    print_array('Computed speed', speed)
    #
    # look up U, V, Speed
    #
    
    wind_u = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['wind_u'].values.copy()
    
    
    wind_v= derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['wind_v'].values.copy()
    
    
    look_speed = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['speed'].values.copy()
    
    print_array('Looked up U', wind_u)
    print_array('Looked up V', wind_v)
    print_array('Looked up speed', look_speed)
    dx0 = dx[1,1]
    dy0 = dy[1,1]
    f0 = f_35deg
    h0 = h500[1,1]
    h1 = h500[1,2]
    h3 = h500[1,0]
    h4 = h500[0,1]
    h2 = h500[2,1]
    laplace0 = h1 + h2 + h3 + h4 - 4*h0
    print("h0 =", h0, "h1=",h1," h2=",h2," h3=", h3," h4=", h4)
    print('laplace0', laplace0)
    d2x0 = d2x[1,1]
    d2y0 = d2y[1,1]
#    zeta1 = (d2x/dx**2 + d2y/dy**2)*GRAVITY/f
    print("d2x0:", d2x0, " d2y0:", d2y0, " dx0:", dx0, " dy0:", dy0)
    zeta1 = (d2x0 / dx0**2 + d2y0 / dy0**2) * GRAVITY /f0
    print('Computed vorticity',zeta1)
