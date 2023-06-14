import pprint
import xarray as xr
import numpy as np
import pandas as pd
from traject_constants import * 
def laplacian(data, dx, dy):
    nrows, ncols = data.shape
    laplace = np.zeros((nrows, ncols))
    print("size: ", nrows, ncols)
    for i in 1, nrows -2:
        for j in 1, ncols -2:
            d2x =  data[i,j+1] - 2*data[i,j] + data[i,j-1]
            d2y = data[i+1,j] - 2*data[i,j] + data[i-1,j]
            laplace[i,j] = d2x/dx**2 + d2y/dy**2

    return laplace

def get_timestamp(start_time, time_step):
    time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_step))
    dt_str = np.datetime_as_string(time_stamp, unit='s')
    return dt_str, time_stamp

def print_array(data):
    print("\t\tLongitude")
    print("lat ", end="")
    for i in range(0,nlon):
        print(f"\t{lon_lin[i]:6.2f}", end="")
    print("")
    for i in range(nlat -1, -1, -1):
            print(lat_lin[i], end="")
            for j in range(0, nlon):
                d = data[i,j]
                print(f"\t{d:6.2f}", end=" ")
            print("")

            


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

dataset = xr.open_dataset("forecasts_baro-spharm_2007-04-15_2007-04-18.nc")
dataset = dataset.roll({"lon": 180})
steps =  dataset.variables['step']
nsteps = len(steps)
step = 7 # 7z on 15th 

x = np.cos(np.deg2rad(lat_lin))*EARTH_RADIUS*np.cos(np.deg2rad(lon_lin))
y = EARTH_RADIUS * np.deg2rad(lat_lin)
dx = np.diff(x)
dy = np.diff(y)
print("dx: ", dx)
print("dy:", dy)
#
timestamps = dataset["time"].values.copy()
forecast_init_time = timestamps[0] # 15Sep 2007
start_time_str, start_time = get_timestamp(forecast_init_time, step)
print("start time: ", start_time_str)
#
# get a window of height data, centered over a jet streak
#
z500 = dataset.sel({
    "time": forecast_init_time,
    "step": steps[step],
    "lat": slice(max_lat, min_lat),
    "lon": slice(180+min_lon, 180+max_lon)
})['z500'].values.copy()

#
# print a table of data
#
#print(z500)
pp = pprint.PrettyPrinter()
array = np.around(z500, decimals=2)
pp.pprint(array)

#print(lat_lin)
#print(lon_lin)
rdiff = np.diff(z500,axis=1)
row_diffs = np.around(np.diff(rdiff, axis=1), decimals=2)
cdiff = np.diff(z500, axis=0)
col_diffs = np.around(np.diff(cdiff, axis=0), decimals=2)
np.set_printoptions(suppress=True)

print(" row diffs:")
pp.pprint(row_diffs)
print(" col diffs:")
pp.pprint(col_diffs)
zeta = laplacian(z500, 1.,1.)
print("laplacian:", zeta.shape)
#pp.pprint(zeta)
print_array(zeta)


