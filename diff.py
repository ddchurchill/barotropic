import numpy as np  
#
# centered differences
#
def centered_diff(z, x, y):
    """
    input;
	z field to differentiate, in both x and y directions
	x matching field of x coordinate
	y matching field of y coordinate
        x and y grids can be variable distances between points.
    output:
	dzdx and dzdy
    
    """
    dzdx = np.zeros_like(z)
    dzdy = np.zeros_like(z)
    nrows, ncols = z.shape
    for i in range(1, ncols -1):
        for j in range(1, nrows -1) :
            dx1 = (x[j,i+1] - x[j,i])
            dzdx1 = (z[j,i+1] - z[j,i])/dx1

            dx2 = (x[j,i] - x[j,i-1])
            dzdx2 = (z[j,i] - z[j,i-1])/dx2
            dzdx[j,i] = (dzdx1 + dzdx2)/2.
            # compute y derivatice
            dy1 = (y[j+1,i] - y[j,i])
            dzdy1 = (z[j+1,i] - z[j,i])/dy1

            dy2 = (y[j,i] - y[j-1,i])
            dzdy2 = (z[j,i] - z[j-1,i])/dy2
            dzdy[j,i] = (dzdy1 + dzdy2)/2.

#
# do edges
#
    for i in range(0,ncols-1):
        dzdx[0,i] = (z[0,i+1] - z[0,i])/(x[0,i+1] - x[0,i])
        dzdx[-1,i] = (z[-1,i+1] - z[-1,i])/(x[-1,i+1] - x[-1,i])
    for j in range(0,nrows -1):
        dzdy[j,0] = (z[j+1,0] - z[j,0])/(y[j+1,0] - y[j,0])
        dzdy[j,-1] = (z[j+1,-1] - z[j,-1])/(y[j+1,-1] - y[j,-1])
#
# now try the numpy way -- not working, shape mismatches.
#
#    dx0 = np.diff(x)
#    dy0 = np.diff(y)
# add padding in first row and column
#    dx =  np.pad(dx0, ((0, 0),(0,1)), mode='constant', constant_values=0)
#    dy =  np.pad(dy0, ((1, 0),(0,0)), mode='constant', constant_values=0)
#    print("centered diff: dx shape:", dx.shape)
#    dzdx[1:-1] = (z[1:] / dx[1:] + z[:-1] / dx[:-1]) / 2
#    dzdx[0] = (z(x[1]) - z(x[0])) / (x[1] - x[0])
#    dzdx[-1] = (z(x[-1]) - z(x[-2])) / (x[-1] - x[-2])
#    dzdy[1:-1] = (z[1:] / dy[1:] + z[:-1] / dy[:-1]) / 2
#    dzdy[0] = (z(y[1]) - z(y[0])) / (y[1] - y[0])
#    dzdy[-1] = (z(y[-1]) - z(y[-2])) / (y[-1] - y[-2])

    return dzdx, dzdy


           

