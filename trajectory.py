#
# Trajectory contains two arrays of variable length, to store the latitude and longitude
# of each point of a trajectory
#
# initialize it with a starting latititude and longitude (in degrees)
class Trajectory:
    def __init__(self, lat0, lon0):
        self.lat = [lat0]
        self.lon = [lon0]
#
# last_lat() returns the last latitude value stored in the object.
# This is referenced by the integration routine euler and huen
#
    def last_lat(self):
        return self.lat[-1]

#
# same for last_lon: returns the last longitude in the array
#
    def last_lon(self):
        return self.lon[-1]
#
# return the array of latitudes
#
    def lats(self):
        return self.lat
    def lons(self):
        return self.lon
    
if __name__ == "__main__":
    a = Trajectory(10., 20.)
    print("a:", a.last_lat(), a.last_lon())
    

    
