from starter1 import *
import yt
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import tools.davetools as dt
import brunt_tools as bt
reload(bt)
import slabimg
reload(slabimg)


L=-0.5; R=0.5; Dx=0.01
x,y,z = np.mgrid[L:R:Dx,L:R:Dx, L:R:Dx]
r2=x**2+y**2+z**2
r = np.sqrt(r2)
thing = np.zeros_like(r)
ok=r>0
thing[ok] = r[ok]**(-1.5)
thing[~ok]==1
thing=np.roll(thing, thing.shape[0]//2, axis=0)
thing=np.roll(thing, thing.shape[1]//2, axis=1)
thing=np.roll(thing, thing.shape[2]//2, axis=2)
thing[:,:,:]=0
thing[5,:,:]=1

Q = np.fft.irfftn(thing[:,:,:thing.shape[2]//2+1])
#fig,axes=plt.subplots(2,2)
#ax0=axes[0][0]

#ax0.imshow(Q.sum(axis=0))
#fig.savefig('plots_to_sort/herp')

ftool=bt.fft_tool(Q)
ftool.do3()
ftool.do2(projax=0)
bt.plot_fft(ftool, outname ='plots_to_sort/cone')
slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
