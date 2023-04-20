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
            dy = (y[j+1,i] - y[j,i])
            dzdx[j,i] = (z[j,i+1] - z[j,i-1])/dx
 
            dzdy[j,i] = (z[j+1,i] - z[j-1,i])/dy

    return dzdx, dzdy



