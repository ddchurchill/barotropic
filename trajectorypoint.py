class TrajectoryPoint:
      # a point on a trajectory
      def __init__(self, lat : float, lon:float, dx:float,
                   dy:float, timestamp):
          self.lat : float = lat
          self.lon : float = lon
          self.dx :float = dx
          self.dy : float = dy
          self.vort : float =0
          self.timestamp = timestamp

