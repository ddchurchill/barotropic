import sys
import pprint
import xarray as xr
import numpy as np
import pandas as pd
from util import get_timestamp
from traject_constants import * 
def getvort_v4(ds, lat, lon, step):

    vort   = ds['abs_vorticity'][step].interp(lat=lat, lon=lon, \
                                              method='slinear',
                                              assume_sorted=True)
    v = vort.values.copy()
#    print("getvort_4: ", lat, lon, step, v)
    return v

def laplacian(data, x_coord, y_coord):
    nrows, ncols = data.shape
    print("Laplacian: shape is ", nrows, ",", ncols)
    laplace = np.zeros((nrows, ncols))
    delsquared = np.zeros((nrows,ncols))
    dx = np.zeros((nrows, ncols))
    ddx = np.zeros((nrows, ncols))
    d2x = np.zeros((nrows, ncols))
    d2y = np.zeros((nrows, ncols))    

    dy = np.zeros((nrows, ncols))
    ddy = np.zeros((nrows, ncols))    
    for i in range(1, nrows -1):
        for j in range( 1, ncols -1):
            dx[i,j] = x_coord[i,j+1] - x_coord[i,j]
            dy[i,j] = y_coord[i+1,j] - y_coord[i,j]
            ddx[i,j] = (data[i,j+1] - data[i,j-1])/(2 * dx[i,j])
            ddy[i,j] = (data[i+1,j] - data[i-1,j])/(2 * dy[i,j])
            d2x[i,j] =  data[i,j+1] - 2*data[i,j] + data[i,j-1]
            d2y[i,j] = data[i+1,j] - 2*data[i,j] + data[i-1,j]
            delsquared[i,j] = d2x[i,j] + d2y[i,j]
            laplace[i,j] = d2x[i,j]/dx[i,j]**2 + d2y[i,j]/dy[i,j]**2            

    return laplace, dx, dy, ddx, ddy, d2x, d2y, delsquared

#
# pretty print an array of data with its title
#
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


if __name__ == "__main__":
    #
    # select a small window of data to print out, and calculate
    # the winds, and vorticity, both computed here and looked up in the
    # dataset file vortdata_steps.nc
    #
    min_lon = -126
    max_lon = -124

    min_lat = 34
    max_lat = 36
    lat0 = 35 # latitude of center of box
    lon0 = -125 # lon of center of box
#    step = 7 # 7z on 15th
    step = 0 # 0z 15th April 2007
    nlat = max_lat - min_lat +1
    
    nlon = max_lon - min_lon +1
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
    print("Dataset steps:", steps)
    
#    step = 7 # 7z on 15th

    
    timestamps = dataset["time"].values.copy()
    forecast_init_time = timestamps[0] # 15Sep 2007
    start_time_str, start_time = get_timestamp(forecast_init_time, step)
    print("Forecast init time: ", forecast_init_time)
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
    
    z500 = np.flip(z500,axis=0) # flip the data as in geowinds4.
    #
    # now get the derived data set
    #
    derived = xr.open_dataset('vortdata_steps.nc')
    
    h500 = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['z500'].values.copy()
    
    h500_a = derived.sel({
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['z500'][step].values.copy()

    looked_rel_vorticity = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['rel_vorticity'].values.copy()
    
    print_array('z500', z500)
    print_array('h500', h500)
    print_array('h500_a', h500_a)
    
    #
    # h500 and z500 should be identical
    #
    mse = np.mean((h500 - z500) ** 2)
    rmse = np.sqrt(mse)
    if rmse != 0:
        print("Error: z500 and h500 differ")
        sys.exit()
    print("Heights match")

    mse = np.mean((h500 - h500_a) ** 2)
    rmse = np.sqrt(mse)
    if rmse != 0:
        print("Error: h500 and h500_a differ")
        sys.exit()
    print("Selected Heights match")
    
    #
    # print a tables of data
    #
    
    
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
    f_lat0 = 2. * OMEGA * np.sin(np.deg2rad(lat0))
    print("Omega: ", OMEGA)
    print("f at center lat: ", f_lat0)
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
    print("Step ", step, " has steps value ", steps[step].values)
    wind_u = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['wind_u'].values.copy()
    
    if wind_u[1,1] == u[1,1] :
        print("U winds match")
    else:
        print("U winds dont match");
        sys.exit()

    wind_v= derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['wind_v'].values.copy()
    
    if wind_v[1,1] == v[1,1] :
        print("V winds match")
    else:
        print("V winds dont match");
        sys.exit()
    
    look_speed = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['speed'].values.copy()
    
    if speed[1,1] == look_speed[1,1]:
        print("Speeds match");
    else:
        print("Speed do not match")
        sys.exit()
        
    print_array('Looked up U', wind_u)
    print_array('Looked up V', wind_v)
    print_array('Looked up speed', look_speed)
    dx0 = dx[1,1]
    dy0 = dy[1,1]
#    f0 = f_35deg
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
    zeta1 = (d2x0 / dx0**2 + d2y0 / dy0**2) * GRAVITY /f_lat0
    print('Computed vorticity',zeta1)
    if zeta1 != looked_rel_vorticity[1,1]:
        print('Computed vorticity and looked up differ')
        sys.exit()
    else:
        print('Computed and looked up vorticity agree')
    #
    # test that the interp function works right
    # 
    #
    zeta2 = getvort_v4(derived, lat0, lon0, step)
    looked_abs_vorticity = derived.sel({
        "step": steps[step],
        "lat": slice(min_lat, max_lat),
        "lon": slice(min_lon, max_lon)
    })['abs_vorticity'].values.copy()

    print("looked up absolute vorticity:", looked_abs_vorticity[1,1])
    print("get_vort finds: ", zeta2)
    if zeta2 == looked_abs_vorticity[1,1]:
        print("getvort_v4 works right")
    else:
        print("getvort is not working right")
        sys.exit()
    # compared with computed absolute vortiticity
    abs_vort_computed = zeta1 + f_lat0
    print("Absolute voriticty computed:", abs_vort_computed)
    if abs_vort_computed == zeta2:
        print("Computed absolute voriticity agreed with getvort_v4")
    else:
        print("Computed absolute voriticitu differs from getvort")
        sys.exit()
    print("All tests passed")
    
        
