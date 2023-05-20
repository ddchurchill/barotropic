from trajectory_v2 import *
from trajectorypoint import TrajectoryPoint 
import matplotlib.pyplot as plt

#
# revised code to compute, plot and analyze trajectories
#
def plot_traject_arrows(trajectory_list):
    """
	plot_traject_arrows - plot the line arrows of a trajectory
	input: list of trajectory of objects
    """	
    #
    # plot trajectories
    #
    # use a different color for nearby trajectories, so they can be easily viewed.
    #
    # the colors of the trajectories cycle through the following list             
    # Can always add more colors to the list as needed.
    colors = ['black', 'red', 'blue', 'green','grey','orange', 'purple']
    color_index = 0
    #
    # repeat for each trajectory
    #
    for trajectory in trajectory_list:
        if trajectory.length == 0 : # skip any empty trajectoriea
            continue
        line_color = colors[color_index % len(colors)]
        color_index += 1
        #
        #	plot each point of the trajectory
        #
        for p in trajectory.points:
    	    plt.arrow(p.lon, p.lat, p.dx, p.dy, \
		length_includes_head=True, head_length=0.8, \
		head_width=0.8, color=line_color)









	

    
