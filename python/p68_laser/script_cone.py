from starter1 import *
import yt
import scipy
import scipy.stats
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import tools.davetools as dt
import brunt_tools as bt
reload(bt)
import slabimg
reload(slabimg)

N = 128
kI = np.fft.fftfreq(N)*N
kx,ky,kz=np.meshgrid(kI,kI,kI)
rrr = np.sqrt(kx**2+ky**2+kz**2)
Ahat = np.zeros_like(rrr*1j)
ok = rrr>0
alpha=-1.5
alphaprime=(alpha-2)/2
Ahat[ok]=rrr[ok]**alphaprime
#Ahat=np.ones([N,N,N])
if 0:
    phi = np.random.random(Ahat.size)
    phi.shape = Ahat.shape
    Ahatmag = np.abs(Ahat)
    Ahat = Ahatmag*np.cos(phi)+Ahatmag*np.sin(phi)*1j
H = Ahat.shape[0]//2
Ahat = Ahat[:,:,:H+1]
Q = np.fft.irfftn(Ahat)
#fig,axes=plt.subplots(2,2)
#ax0=axes[0][0]

#ax0.imshow(Q.sum(axis=0))
#fig.savefig('plots_to_sort/herp')

if 1:
    #FFT = np.fft.fftn(Q)
    #power = FFT*np.conj(FFT)
    power = (Ahat*Ahat.conj()).real
    kbins = np.sqrt( np.unique(kI**2))
    kbins.sort()
    kcen = 0.5*(kbins[1:]+kbins[:-1])

    kvals = rrr[:,:,:H+1]
    TheY, bin_edge, count = scipy.stats.binned_statistic(kvals.flatten(), power.flatten(), bins=kbins, statistic='sum')

    TheX = kcen
    ok = (TheX>0)*(TheY>0)
    pfit = np.polyfit( np.log(TheX[ok]), np.log(TheY[ok]), 1)
    print('I get',pfit)

if 1:
    bt.plot_fft(ftool, outname ='plots_to_sort/cone')
    slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
    fig,ax=plt.subplots(1,1)
    Q=ftool.power2.real
    Q=np.roll(Q,Q.shape[0]//2,axis=0)
    Q=np.roll(Q,Q.shape[1]//2,axis=1)
    ki=np.fft.fftfreq(Q.shape[0])
    kx,ky=np.meshgrid(ki,ki)
    r2d = np.sqrt(kx**2+ky**2)
    norm=mpl.colors.LogNorm(vmin=Q[Q>0].min(),vmax=Q.max())
    denom = np.ones_like(r2d)
    ok = r2d>0
    denom[ok] = r2d[ok]**alpha
    ax.imshow(Q,norm=norm,cmap='gray')
    fig.savefig('plots_to_sort/power2d')

if 'ftool' not in dir():
    ftool=bt.fft_tool(Q)
    ftool.do3()
    ftool.do2(projax=0)
if 1:
    slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
    if 'slab' not in dir():
        slab=slabimg.slab(ftool, projax=0)
    slabimg.number_checker(ftool, outname = "plots_to_sort/cone_number")
    slabimg.plot(ftool,outname='plots_to_sort/cone_diss')
