import pprint
import xarray as xr
import numpy as np
import pandas as pd
from traject_constants import * 
def laplacian(data, x, y):
    nrows, ncols = data.shape
    laplace = np.zeros((nrows, ncols))
    dx = np.zeros((nrows, ncols))
    ddx = np.zeros((nrows, ncols))
    dy = np.zeros((nrows, ncols))
    ddy = np.zeros((nrows, ncols))    
    print("size: ", nrows, ncols)
    for i in 1, nrows -2:
        for j in 1, ncols -2:
            dx[i,j] = x[i,j+1] - x[i,j]
            dy[i,j] = y[i+1,j] - y[i,j]
            ddx[i,j] = (data[i,j+1] - data[i,j-1])/(2 * dx[i,j])
            ddy[i,j] = (data[i+1,j] - data[i-1,j])/(2 * dy[i,j])
            d2x =  data[i,j+1] - 2*data[i,j] + data[i,j-1]
            d2y = data[i+1,j] - 2*data[i,j] + data[i-1,j]
            laplace[i,j] = d2x/dx[i,j]**2 + d2y/dy[i,j]**2

    return laplace, dx, dy, ddx, ddy

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
dataset = dataset.roll({"lon": 180})

steps =  dataset.variables['step']
nsteps = len(steps)
step = 7 # 7z on 15th 

timestamps = dataset["time"].values.copy()
forecast_init_time = timestamps[0] # 15Sep 2007
start_time_str, start_time = get_timestamp(forecast_init_time, step)
print("start time: ", start_time_str)
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

rel_vorticity = derived.sel({
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


zeta, dx, dy, ddx, ddy  = laplacian(h500, x, y)

print_array("Dx", dx)
print_array('Dy', dy)
print_array('ddx', ddx)
print_array('ddy', ddy)
print_array('Laplacian', zeta)
print_array('coriolis', f)

# compute vorticity, multiply by gravity, divide by f
vorticity = zeta * GRAVITY / f
print_array('vorticity', vorticity)
f_35deg = 2. * OMEGA * np.sin(np.deg2rad(35.))
print("Omega: ", OMEGA)
print("f 35 deg: ", f_35deg)
print_array('Rel vorticity', rel_vorticity)

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
