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
        
    def plot_4panel(self, dataset, time_index, time_str, show_plot):

        zeta = dataset['rel_vorticity'][time_index].values.copy()
        abs_vor = dataset['abs_vorticity'][time_index].values.copy()
        z500 = dataset['z500'][time_index].values.copy()
        speed = dataset['speed'][time_index].values.copy()

        u = dataset['wind_u'][time_index].values.copy()
        v = dataset['wind_v'][time_index].values.copy()
    
        fig, ax = plt.subplots(2,2,\
                               figsize=(20, 10), \
                               subplot_kw={'projection': ccrs.PlateCarree()})
        fig.suptitle('Winds, Heights, Vorticity at ' + time_str)

        # Draw parallels and meridians
        parallels = range(min_lat, max_lat, 10)
        meridians = range(min_lon, max_lon, 10)

        for axis in ax.flat:
            axis.set_extent([min_lon, max_lon, min_lat, max_lat])

        z500_con = ax[0,0].contourf(lon_lin, lat_lin,z500,\
                            vmin=5000., vmax=6000., \
                            levels=16, cmap='jet', \
                            transform=ccrs.PlateCarree())

        z500_bar = plt.colorbar(z500_con, ax=ax[0,0], orientation='horizontal')
        units  = r'$(m)$'
        z500_bar.set_label('500 mb heights ' + units)
        z500_con.set_clim(vmin=5000, vmax=6000)        
        ax[0,0].set_title('Winds and Heights ' )

        # draw the winds
        m = Basemap(projection='cyl', llcrnrlat=min_lat, \
                urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon)

        x, y = m(lon, lat)

        ax[0,0].quiver(x[::5, ::5], y[::5, ::5], u[::5, ::5], v[::5, ::5], \
                       scale=3000, color='black')

        # plot relative vorticity
        max_vort = 5.e-4  # uppeer limit for vorticity plot
        zeta_con = ax[0,1].contourf(lon_lin, lat_lin, zeta, cmap='jet', \
                                    levels = 16)
# vmin=-max_vort, vmax=max_vort, extend='neither')

        zeta_bar = plt.colorbar(zeta_con, ax=ax[0,1], orientation='horizontal')
        ax[0,1].set_title('Relative Vorticity ' )
        units  = r'$(s^{-1})$'
        zeta_bar.set_label('Vorticity ' + units)
#        zeta_con.set_clim(vmin=-1.e-4, vmax=2.e-4)                 
        # plot absolute vorticity

        abs_vor_con = ax[1,1].contourf(lon_lin,lat_lin,abs_vor, \
                                       levels=16, cmap='jet',
                                vmin=0., vmax=max_vort, extend='neither')

        abs_vor_bar = plt.colorbar(abs_vor_con, ax=ax[1,1], orientation='horizontal')
        ax[1,1].set_title('Absolute Vorticity ')
        units  = r'$(s^{-1})$'
        abs_vor_bar.set_label('Vorticity ' + units)
        abs_vor_con.set_clim(vmin=0., vmax=max_vort)        

        # plot wind speed
        speed_con = ax[1,0].contourf(lon_lin,lat_lin,speed, \
                                     vmin=0, vmax=80.,
                                     levels=16, cmap='jet')
        ax[1,0].set_title('Wind Speed ' )
        speed_bar = plt.colorbar(speed_con, ax=ax[1,0], orientation='horizontal')
        units  = r'$(m \,s^{-1})$'
        speed_bar.set_label('Speed ' + units)
        speed_con.set_clim(vmin=0, vmax=80)        
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
            axis.set_ylim(bottom=min_lat, top=max_lat)
        plt.tight_layout() # this is needed to prevent overlapping figures.

        filename = 'wind-vort-' + time_str + '.png'
        plt.savefig(filename)

        if show_plot:
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
    p = wxplots() # class containing plotting routine
    show_plot = False
    for time_index in range(0,48):
        time_stamp = np.datetime64(start_time + pd.Timedelta(hours=time_index))
        time_str = np.datetime_as_string(time_stamp, unit='s')
        print("plotting " , time_str)
        p.plot_4panel(data, time_index, time_str, show_plot)

