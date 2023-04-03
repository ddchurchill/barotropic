import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#
# perform one step of Huen's integration
#

gravity = 9.81 # m/s/s
earth_radius = 6371.e3 # meters
# Create a grid of latitude and longitude values
# 1 degree grids
grid_spacing = 1. # 1 degree grid spacing
nlat, nlon = 140, 360
lat_lin = np.linspace(-70, 70, nlat)
mu = np.sin(lat_lin * np.pi/180.) # used in equation for geostrophic east-west wind
dmu = np.gradient(mu) 
lon_lin = np.linspace(-180, 180, nlon)
lon, lat = np.meshgrid(lon_lin, lat_lin)
phi_rad = np.array(lat_lin * np.pi/180.)  # latitude in radians
lambda_rad = np.array(lon_lin * np.pi/180.)  # longitude in radians

# Define the geopotential field with two ridges in the northern hemisphere and two ridges in the southern hemisphere
# north pacific
decay = 10
r1 = np.exp(-((lat - 30) / decay) ** 2 - ((lon -170) / decay) ** 2)

r2 = np.exp(-((lat + 30) / decay) ** 2 - ((lon + 120) / decay) ** 2)
# north atlantic low
r3 = - np.exp(-((lat - 30) / decay) ** 2 - ((lon + 40) / decay) ** 2)
r4 = np.exp(-((lat + 30) / decay) ** 2 - ((lon + 10) / decay) ** 2)

geopot = ((r1 + r2 + r3 + r4) * 100 + 5000) * gravity


# Compute the geostrophic winds
# This works
#dzdx = np.gradient(geopot, axis=1) / np.gradient(lon * np.pi / 180, axis=1)
#
# try dzdz as second argument. Does not work.
#
lon_grad = np.gradient(lon_lin* np.pi/180) # 1-D array longitudes
#dzdx = np.gradient(geopot, lon_grad, axis=1)

deltay = np.pi/180./earth_radius
deltax = np.pi/180./earth_radius # wrong - needs cos phi
lat_grad = np.gradient(lat * np.pi / 180)
dzdy = np.gradient(geopot, axis=0) / np.gradient(lat * np.pi / 180, axis=0)
#
#
# calculate dzdy using constant difference argument of 1 deg resolution
#
dzdy = np.gradient(geopot,grid_spacing * np.pi / 180, axis=0) 

f = 2 * np.pi / 86400 * np.sin(lat * np.pi / 180)

#
#
# dzdx is change of phi with respect to lambda
#
dzdx = np.gradient(geopot, grid_spacing * np.pi/180, axis=1)

#dzdmu = np.gradient(geopot, dmu, axis=0)
ug = - dzdy / (f* earth_radius)
vg = dzdx/ (f * earth_radius * np.cos(lat * np.pi/180)) 

# compute geostrophic wind speed
geostrophic_wind_speed = np.sqrt(ug**2 + vg**2)
#
data_range = np.ptp(geostrophic_wind_speed)
print("Range:", data_range)
mean = np.mean(geostrophic_wind_speed)
print("Mean speed:", mean)
maxspeed = np.max(geostrophic_wind_speed)
print("max speed:", maxspeed)

# Create a new figure
fig = plt.figure(figsize=(12, 6))

# Create a new map projection
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)

# Draw the continents and coastlines
m.drawcoastlines(linewidth=0.5,color='white')
m.drawcountries(linewidth=0.5, color='white')

# Draw the geopotential field
x, y = m(lon, lat)
m.contourf(x, y, geopot, cmap='jet', levels=50)

# Draw the geostrophic wind vectors
m.quiver(x[::5, ::5], y[::5, ::5], ug[::5, ::5], vg[::5, ::5], scale=1000, color='white')

# Add a colorbar and title
plt.colorbar(label='Geopotential')
plt.title('Geopotential and Geostrophic Wind Vectors')
plt.show()



f = 2 * np.pi / 86400 * np.sin(lat * np.pi / 180)

# Create a separate plot for the geostrophic wind speed
fig2 = plt.figure(figsize=(12, 8))

# Create a new map projection
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)

# Draw the continents and coastlines in white
m.drawcoastlines(linewidth=0.5, color='white')
m.drawcountries(linewidth=0.5, color='white')

m.contourf(x,y,geostrophic_wind_speed, cmap='jet',levels=50)
# Draw the geostrophic# Show the plot
# Add a colorbar and title
plt.colorbar(label='Wind speed')
plt.title('Geostrophic Wind Speed')
plt.show()

#
# plot the absolute voriticity
#

deltax = 1/(earth_radius * np.cos(lat))
deltay = 1./earth_radius
#d2zdx2 = np.gradient(dzdx,axis=1)/np.gradient(lon * np.pi/180, axis = 1)
#d2zdy2 = np.gradient(dzdy,axis=0)/np.gradient(lat * np.pi/180, axis = 0)

coslat = np.cos(lat * np.pi/180)
#vorticity = (d2zdx2 + d2zdy2) # + f
#dvdy = np.gradient(vg, axis=0)/np.gradient(lat * np.pi/180,axis=0)/earth_radius
#dudx = np.gradient(ug,axis=1)/np.gradient(lon * np.pi/180, axis=1) / (earth_radius*coslat) 
dvdx = np.gradient(vg, axis=0)/np.gradient(lon * np.pi/180, axis=1) / (earth_radius*coslat)
dudy = np.gradient(ug,axis=1)/np.gradient(lon * np.pi/180, axis=1) / (earth_radius*coslat) 

vorticity = dvdx - dudy  + f
#
# use laplacian to get vorticity
#

dphi = np.gradient(phi_rad)
dlambda = np.gradient(lambda_rad)
dzdphi = np.gradient(geopot,axis=0)/dphi
dzdlambda = np.gradient(geopot,axis=1)/dlambda
d2zdphi2 = np.gradient(dzdphi, axis=0)/dphi
d2zdlambda2 = np.gradient(dzdlambda, axis=1)/dlambda

vorticity = d2zdphi2/earth_radius**2 + d2zdlambda2/(earth_radius*np.cos(phi_rad))**2

#vorticity = np.gradient(np.gradient(geopot,axis=0),axis=0)/(earth_radius **2) \
#    np.gradient(np.gradient(geopot,axis=1),axis=1)/(earth_radius * cos_lat)**2 

print("max vort: ", np.max(vorticity))
print("min vort: ", np.min(vorticity))
      
fig3 = plt.figure(figsize=(12,8))

m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)

# Draw the continents and coastlines in white                                                                            
m.drawcoastlines(linewidth=0.5, color='white')
m.drawcountries(linewidth=0.5, color='white')

m.contourf(x,y,vorticity, cmap='jet',levels=50)
# Add a colorbar and title                                                                                               
plt.colorbar(label='vorticity')
plt.title('absolute vorticity')
plt.show()
#
# now calculate trajectories
# start with Euler's method
#
npart = 10 # number of parcel trajectories
nt = 10 # number of time steps
#
# set initial lat and lon of trajector
#
init_lat = 35.
init_lon = -40.
# Trajectory calculation using Euler's method
lat_traj = np.zeros((npart, nt))  # Array to store latitude trajectory
lon_traj = np.zeros((npart, nt))  # Array to store longitude trajectory
lat_traj[:, 0] = init_lat
lon_traj[:, 0] = init_lon
coriolis = 2 * np.pi / 86400
#
# repeat for each time step
#
for t in range(1, nt):
    # Compute Coriolis parameter at current latitude
    f_lat = coriolis * np.sin(np.radians(lat_traj[:, t-1]))
    
    # Compute velocity at previous location
    u_prev = np.interp(lon_traj[:, t-1], lon, ug[np.argmin(np.abs(lat_traj[:, t-1] - lat))])
    v_prev = np.interp(lon_traj[:, t-1], lon, vg[np.argmin(np.abs(lat_traj[:, t-1] - lat))])
    
    # Compute velocity at current location
    u = u_prev
    v = v_prev + f_lat * u_prev * dt
    
    # Compute new position using Euler's method
    lat_traj[:, t] = lat_traj[:, t-1] + v * dt
    lon_traj[:, t] = lon_traj[:, t-1] + u * dt

# Plotting the trajectories
plt.figure(figsize=(8, 6))
plt.title('Particle Trajectories')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
for i in range(npart):
    plt.plot(lon_traj[i], lat_traj[i], linewidth=2)
plt.grid()
plt.show()
