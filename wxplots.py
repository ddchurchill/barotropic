import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy
from mpl_toolkits.basemap import Basemap

from constants import *
class wxplots:
    def __init__(self):
        self.a = 0
        
    def plot_4panel(self, dataset, time_index, time_str):

        zeta = dataset['rel_vorticity'][time_index].values.copy()
        abs_vor = dataset['abs_vorticity'][time_index].values.copy()
        z500 = dataset['z500'][time_index].values.copy()
        speed = dataset['speed'][time_index].values.copy()
        u = dataset['wind_u'][time_index].values.copy()
        v = dataset['wind_v'][time_index].values.copy()
        x = np.cos(lat_rad)*EARTH_RADIUS*lon_rad
        y = EARTH_RADIUS * lat_rad

    
        fig, ax = plt.subplots(2,2,\
                               figsize=(20, 10), \
                               subplot_kw={'projection': ccrs.PlateCarree()})

        # Draw parallels and meridians
        parallels = range(min_lat, max_lat, 10)
        meridians = range(min_lon, max_lon, 10)

        for axis in ax.flat:
            axis.set_extent([min_lon, max_lon, min_lat, max_lat])

        con1 = ax[0,0].contourf(lon_lin, lat_lin,z500,\
                            vmin=5000., vmax=6000., \
                            levels=16, cmap='jet', \
                            transform=ccrs.PlateCarree())

        cbar = plt.colorbar(con1, ax=ax[0,0], orientation='horizontal')
        cbar.set_label('(m)')
        
        ax[0,0].set_title('Winds and Heights ' + time_str)

        # draw the winds
        m = Basemap(projection='cyl', llcrnrlat=min_lat, \
                urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon)

        x, y = m(lon, lat)

        ax[0,0].quiver(x[::5, ::5], y[::5, ::5], u[::5, ::5], v[::5, ::5], \
                       scale=3000, color='black')


        con2 = ax[0,1].contourf(lon_lin, lat_lin, zeta, cmap='jet', \
                                 levels = 16)
#                                 vmin=0., vmax=3.e-4, extend='neither')

        plt.colorbar(con2, ax=ax[0,1], orientation='horizontal')
        ax[0,1].set_title('Relative Vorticity ' + time_str)

         
        con3 = ax[1,1].contourf(lon_lin,lat_lin,abs_vor, \
                                levels=16, cmap='jet', \
                                vmin=0., vmax=3.e-4, extend='neither')

        plt.colorbar(con3, ax=ax[1,1], orientation='horizontal')
        ax[1,1].set_title('Absolute Vorticity ' + time_str)
    



        con4 = ax[1,0].contourf(lon_lin,lat_lin,speed, levels=16, cmap='jet')
        ax[1,0].set_title('Wind Speed ' + time_str)
        plt.colorbar(con4, ax=ax[1,0], orientation='horizontal')

        # Draw the continents and coastlines in white
        for axis in ax.flat:
            axis.coastlines(linewidth=0.5, color='black')
            axis.add_feature(cartopy.feature.BORDERS, linewidth=0.5, \
                             edgecolor='black')
            axis.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5,
                           linestyle='--')
            axis.set_xticks(meridians, crs=ccrs.PlateCarree())
            axis.set_yticks(parallels, crs=ccrs.PlateCarree())
            axis.xaxis.set_ticklabels([])
            axis.yaxis.set_ticklabels([])

        plt.tight_layout() # this is needed to prevent overlapping figures.

            #    plt.savefig(filename)

        plt.show()
        plt.close()
# main program
if __name__ == "__main__":
    dataset_file = 'vortdata_steps.nc'

    data = xr.open_dataset(dataset_file)
    print("Data read in from ", dataset_file)
    time_index = 1
    start_time = data["time"].values
    print("Start time: ", start_time)
    time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_index))
    time_str = np.datetime_as_string(time_stamp, unit='s')
    p = wxplots()
    p.plot_4panel(data, time_index, time_str)

