from starter1 import *
import yt
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import tools.davetools as dt
import brunt_tools as bt
reload(bt)
import slabimg
reload(slabimg)

N = 128
L=-0.5; R=0.5; Dx=1/N
kI = np.fft.fftfreq(N)*N
kx,ky,kz=np.meshgrid(kI,kI,kI)
r2 = kx**2+ky**2+kz**2
Ahat = np.zeros_like(r2*1j)
ok = r2>0
Ahat[ok]=r2[ok]**-1.5
if 0:
    phi = np.random.random(Ahat.size)
    phi.shape = Ahat.shape
    Ahatmag = np.abs(Ahat)
    Ahat = Ahatmag*np.cos(phi)+Ahatmag*np.sin(phi)*1j
H = Ahat.shape[0]//2
Ahat = Ahat[:,:H+1]
Q = np.fft.irfftn(Ahat)
#fig,axes=plt.subplots(2,2)
#ax0=axes[0][0]

#ax0.imshow(Q.sum(axis=0))
#fig.savefig('plots_to_sort/herp')

ftool=bt.fft_tool(Q)
ftool.do3()
ftool.do2(projax=0)
bt.plot_fft(ftool, outname ='plots_to_sort/cone')
slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
if 'slab' not in dir():
    slab=slabimg.slab(ftool, projax=0)
slabimg.number_checker(ftool, outname = "plots_to_sort/cone_number")
