class TrajectoryPoint:
      # a point on a trajectory
      def __init__(self, lat, lon, dx, dy, timestamp):
          self.lat = lat
          self.lon = lon
          self.dx = dx
          self.dy = dy
          self.vort =0
          self.timestamp = timestamp

