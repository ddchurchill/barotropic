import numpy as np
import xarray as xr
from huen_v4 import getvort_v4 # looks up vorticity values
#
# read in the trajectory file, and print out the values
# Check that the vorticity stored matches what is in the file
# vortdata_steps.nc
#
trajectory_file = 'trajectories.npz' # numpy compressed file

npdata = np.load(trajectory_file,allow_pickle=True)
trajectories = npdata['trajectories']
print("Read in trajectory data from ", trajectory_file)
dataset = xr.open_dataset('vortdata_steps.nc')

failed = False
for i, t in enumerate(trajectories):
    print("Trajectory ", i," start:", t.start_time, " stop:", t.stop_time)
    for j, p in enumerate(t.points):
        if np.isnan(p.vort):
            continue  # skip NaN values
        step = j # time step is just an iterator
        v = getvort_v4(dataset, p.lat, p.lon, step)
        error = np.abs(v - p.vort)
        if error > 1.e-6:
            print("Point ", j,"time:", p.timestamp,"vort:", p.vort)
            print("\t\t\tlat:", p.lat, "\t lon:", p.lon  )
            print("\t\t\tReference vorticity: ", v)
            failed = True
if failed:
    print("Test failed")
else:
    print("Test passes")
    
        
