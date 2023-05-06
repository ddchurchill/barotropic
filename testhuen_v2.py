min_lon = -60
max_lon = -40
min_lat = 10
max_lat = 30
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from velocity import Velocity as vel
from trajectory_v2 import Trajectory_v2
from huen_v2 import *
import traject_constants
#
    
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
        return vel(-10,10)
    if(time_step == 7):
        return vel(-20, 0)
    if(time_step == 8):
        return vel(-20, 0)
    if(time_step == 9):
        return vel(-10, -10)
    if(time_step == 10):
        return vel(-10, -10)
    if(time_step == 11):
        return vel(0, -20)
    if(time_step == 12):
        return vel(0, -20)
    if(time_step == 13):
        return vel(10,-10)
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
    if(time_step == 11):
        return vel(0, 20)
    if(time_step == 12):
        return vel(-10, 0)
    if(time_step == 13):
        return vel(-10, 0)
    return vel(0,20)

        

    
#
# nowind returns null values for winds, as when lat and lon are out of bounds
#
def  nowind(lat, lon, time_step):
    return None

def printtraj(trajectory, label):
    print(label)
    fig4 = plt.figure(figsize=(12,8))


    # Draw the continents and coastlines in white                            
    m.drawparallels(range(min_lat,max_lat, 10), labels=[1,0,0,0])
    m.drawmeridians(range(min_lon, max_lon, 10), linewidth=1, labels=[0,0,0,1])
    for p in trajectory.points:


        print(f"{p.lat:.2f}", f"{p.lon:.2f}", \
        f"{p.dx:.2f}", f"{p.dy:.2f}")
        plt.arrow(p.lon, p.lat, p.dx, p.dy, color='red',\
              length_includes_head=True, head_length=0.1, head_width=0.5)
    plt.title(label)
    plt.show()

def main():
    nt = 13
    deltat = 3600  # 1 hour integration period
    lat0 = 20. # degrees
    lon0 = -50. # degrees starting point

##    traj = huen_v2(anticyclone, lat0, lon0, deltat, nt)



    np.set_printoptions(precision=6)
    model = huen_v2

    traj = model(westwind, lat0, lon0, deltat, nt)
    printtraj(traj,"westwind")

    traj = model(northwind, lat0, lon0, deltat, nt)
    printtraj(traj, "northwind")

    traj = model(cyclone, lat0, lon0, deltat, nt)
    printtraj(traj, "cyclone")

    traj = model(anticyclone, lat0, lon0, deltat, nt)
    printtraj(traj, "anticyclone")
    

#
# now test when no winds are found
#
#    lat_traj, lon_traj = huen(nowind, lat0, lon0, deltat, nt)
#    print("no wind lat: ",lat_traj)
#    print("no wind lon: ", lon_traj)
    
if __name__ == "__main__":
    m = Basemap(projection='cyl', llcrnrlat=min_lat, urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon)
    
    main()          
