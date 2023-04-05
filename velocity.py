import numpy as np
class Velocity:
    def __init__(self, u, v):
        self.u = u
        self.v = v

    def speed(self):
        s = (self.u**2 + self.v**2)**0.5  # calculate the speed
        return s

    def direction(self): # returns direction in degrees. (Not from North)
        d = np.arctan2(self.v, self.u)
        return np.degrees(d)
        
    def __str__(self):
        return f"Velocity: ({self.u}, {self.v}), speed: {self.speed():.2f}"

if __name__ == "__main__":
        v = Velocity(10,10)
        print("Speed: ", v.speed())
        print("Direction:", v.direction())
        
