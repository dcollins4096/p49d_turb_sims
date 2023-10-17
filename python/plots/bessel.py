
from GL import *
import scipy

x,y = np.mgrid[0:10:0.1,0:10:0.1]
a = scipy.special.j0(x*y)
fig,axes=plt.subplots(1,2)
ax=axes[0];ax1=axes[1]
#plot=ax.imshow(a)
plot=ax.pcolormesh(x,y,a, shading='nearest')

z = np.linspace(0,10,100)
ax1.plot( z, scipy.special.j0(z))
fig.colorbar(plot,ax=ax)
fig.savefig('plots_to_sort/bessel')
plt.close('all')

fig,ax0=plt.subplots(1,1)
k = np.linspace(1,256,256)
k/=k.max()
y1=k**(-1.5)
y3=k**(-3)
y2 = y1*scipy.special.j0(10*k)
ax0.plot(k, y1,c='r')
ax0.plot(k,y2,c='g')
ax0.plot(k,y3,c='b')
ax0.set(xscale='log',yscale='log')
fig.savefig('plots_to_sort/powers')

