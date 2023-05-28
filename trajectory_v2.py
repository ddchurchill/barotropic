class Trajectory_v2:  # revised trajectory class
    def __init__(self):
        self.points = [TrajectoryPoint() for _ in range(0)]
        self.length = 0
    def __init__(self, lat0, lon0, deltat, start_time, stop_time):
        self.points = [TrajectoryPoint() for _ in range(0)]
        self.length = 0
        self.start_lat = lat0
        self.start_lon = lon0
        self.time_step = deltat # in seconds
        self.last_lat = lat0 # ending latitude of trajectory
        self.last_lon = lon0 # ending longitude of trajectory
        self.start_time = start_time # starting time stamp of trajectory
        self.stop_time = stop_time # ending time stamp of trajectory
        


      	  
