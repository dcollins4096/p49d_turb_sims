import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def do_stuff(xdata, ydata, p0):
    def fitFunc(x, a, b, c, d):
    	return a+ b*x[0] + c*x[1] + d*x[0]*x[1]


    fitParams, fitCovariances = curve_fit(fitFunc, xdata, ydata, p0)
    print(' fit coefficients:\n', fitParams)
    return fitParams, fitCovariances

x_3d = np.array([[1,2,3,4,6],[4,5,6,7,8]])
y = x_3d[1,:]
p0 = [5.11, 3.9, 5.3, 2]
fitParams, fitCovariances= do_stuff(x_3d,y, p0)

x=np.row_stack([the_x.flatten(),the_y.flatten()])
y=matrix.flatten()
p0 = [0,1,1,0]
fitParams, fitCovariances= do_stuff(x,y, p0)

fig, axx=plt.subplots(2,1)
ax=axx[0]
ax1=axx[1]
for nsub in range(3):
    ax.plot( the_x[:,nsub], matrix[:,nsub], c='rgb'[nsub],marker='*')
    f = fitParams
    amps = f[0] + f[1]*the_x[:,nsub] + f[2]*the_y[:,nsub]+f[3]*the_x[:,nsub]*the_y[:,nsub]
    ax.scatter( the_x[:,nsub], amps,c='k')
for nsub in range(4):
    ax1.plot( the_y[nsub,:], matrix[nsub,:])
    amps = f[0] + f[1]*the_x[nsub,:] + f[2]*the_y[nsub,:]+f[3]*the_x[nsub,:]*the_y[nsub,:]
    ax1.scatter( the_y[nsub,:] , amps)
fig.savefig('%s/bilinear.png'%plotdir)
