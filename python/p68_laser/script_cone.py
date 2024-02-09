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
alphaT=-1.5
alphaV=alphaT-2
kmin=2
kmax=-2
Q = bt.fake_powerlaw(N,alphaT,kmin,kmax)
#fig,axes=plt.subplots(2,2)
#ax0=axes[0][0]

#ax0.imshow(Q.sum(axis=0))
#fig.savefig('plots_to_sort/herp')



if 'ftool' not in dir():
    ftool=bt.fft_tool(Q)
    ftool.do3()
    ftool.do2(projax=0)

    #if 'slab' not in dir():
    #    slab=slabimg.slab(ftool, projax=0)

if 1:
    bt.plot_fft(ftool, outname ='plots_to_sort/cone')
    slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
    slabimg.number_checker(ftool, outname = "plots_to_sort/cone_number")

if 1:
    FFT = np.fft.fftn(Q)
    power = (FFT*np.conj(FFT)).real
    #power = (Ahat*Ahat.conj()).real
    kI = np.fft.fftfreq(N)
    kx,ky,kz=np.meshgrid(kI,kI,kI)
    rrr = np.sqrt(kx**2+ky**2+kz**2)
    H = Q.shape[0]//2
    kbins = np.sqrt( np.unique(kI**2))
    rmin=kbins[kmin]
    rmax=kbins[kmax]
    
    kbins.sort()

    kcen = 0.5*(kbins[1:]+kbins[:-1])

    kvals = rrr[:,:,:]
    TheY, bin_edge, count = scipy.stats.binned_statistic(kvals.flatten(), power.flatten(), bins=kbins, statistic='sum')

    TheX = kcen
    ok = (TheX>0)*(TheY>0)
    pfit = np.polyfit( np.log(TheX[ok]), np.log(TheY[ok]), 1)
    print('I get',pfit)

if 1:
    Nzones = ftool.power.size
    dk = bin_edge[1:]-bin_edge[:-1]
    I = (ftool.power_1d3).sum()
    Ialso = 4*np.pi/(alphaV+3)*(rmax**(alphaV+3)-rmin**(alphaV+3))
    I3 = ftool.power.real.sum()
    I2d1 = 2*np.pi/(alphaV+2)*(rmax**(alphaV+2)-rmin**(alphaV+2))
    I2d2 = (ftool.power_1d2.real).sum()
    I2d3 = (ftool.power[0,:,:].real).sum()*128
    print("Spect I %.3f"%I.real)
    print("Analyt  %.3f"%Ialso)
    print("Cube    %.3f"%(I3))
    print("---")
    print("2d Ann  %.3f"%(I2d1))
    print("2d Coun %.3f"%(I2d2))
    print("2d 3d   %.3f"%(I2d3))
    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]
    ax0.imshow(power.sum(axis=0))
    ax1.imshow(ftool.power.real.sum(axis=0))
    if 0:
        ok = (rrr>0)*(q>0)
        q = ftool.power.real*rrr**2
        ax2.scatter(rrr[ok].flatten(),q[ok].flatten())
        ax2.set(yscale='log')
        q = ftool.power.real
        ax3.scatter(rrr[ok].flatten(),q[ok].flatten())
        #ax3.set(yscale='log')
    fig.savefig('plots_to_sort/cone_power')


if 0:
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

if 0:
    slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
    if 'slab' not in dir():
        slab=slabimg.slab(ftool, projax=0)
    slabimg.number_checker(ftool, outname = "plots_to_sort/cone_number")
    slabimg.plot(ftool,outname='plots_to_sort/cone_diss')
