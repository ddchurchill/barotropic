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
    output:
	dzdx and dzdy
    
    """
    dzdx = np.zeros_like(z)
    dzdy = np.zeros_like(z)
    nrows, ncols = z.shape
    for i in range(2, ncols -1):
        for j in range(2, nrows -1) :
            dx = (x[j,i+1] - x[j,i-1])

            dy = (y[j+1,i] - y[j-1,i])
            dzdx[j,i] = (z[j,i+1] - z[j,i-1])/dx
 
            dzdy[j,i] = (z[j+1,i] - z[j-1,i])/dy

    return dzdx, dzdy
#
# second order centered difference
#
def centered_diff2(z,x,y):
    print("Centered_diff2")
    d2zdx2 = np.zeros_like(z)
    d2zdy2 = np.zeros_like(z)
    nrows, ncols = z.shape
    for i in range(2, ncols -1):
        for j in range(2, nrows -1) :
            dx = (x[j,i+1] - x[j,i-1]) / 2.

            dy = (y[j+1,i] -  y[j-1,i]) / 2.
            d2zdx2[j,i] = (z[j,i+1] -2*z[j,i] + z[j,i-1])/dx**2
 
            d2zdy2[j,i] = (z[j+1,i] -2* z[j][i] + z[j-1,i])/dy**2

    return d2zdx2, d2zdy2


           

