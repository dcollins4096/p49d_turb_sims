
from starter1 import *
import yt
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
import tools.pcolormesh_helper as pch
import tools.davetools as dt
reload(pch)
reload(bt)
plt.close('all')

sim='6_2'
frame=51
sim='half_half'
frame=31

#get cubes; rho_full is straight off disk.  rho is downsampled by 2
if 'rho1' not in dir():
    print('load and cg')
    rho_full, rho_256, rho1 = bt.get_cubes(sim,frame, do_rho_4=True)
    print(sim,frame)
    ft1=bt.fft_tool(rho1)
    ft1.do3()
    ft1.do2(projax=0)

if 'rho2' not in dir() or True:
    N = 128
    alphaT=-1.5
    alphaV=alphaT-2
    kslice=slice(1,-2)
    kslice=slice(None)
    kslice=None
    rho2 = bt.fake_powerlaw(N,alphaT,kslice, Amplitude=1e3, phase=True, rando=True)

    ft2=bt.fft_tool(rho2)
    ft2.do3()
    ft2.do2(projax=0)

if 1:
    kI = np.fft.fftfreq(N)
    kx,ky,kz=np.meshgrid(kI,kI,kI)
    rrr = np.sqrt(kx**2+ky**2+kz**2)
    theta = np.zeros_like(rrr)
    ok = rrr>0
    theta[ok] = np.arccos( kz[ok]/rrr[ok])
    phi = np.arctan2(kx,ky)

if 0:
    #Plot: phases
    fig,axes=plt.subplots(4,2,figsize=(12,12))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ax4=axes[2][0];ax5=axes[2][1]
    ax6=axes[3][0];ax7=axes[3][1]

    N = ft1.fft3.shape[0]
    FFT1=ft1.fft3[N//2,:,:]
    FFT2=ft2.fft3[N//2,:,:]
    p1=ax0.imshow(FFT1.real)
    fig.colorbar(p1,ax=ax0)
    p1=ax1.imshow( FFT2.real)
    fig.colorbar(p1,ax=ax1)
    p1=ax2.imshow( FFT1.imag)
    fig.colorbar(p1,ax=ax2)
    p1=ax3.imshow( FFT2.imag)
    bins=np.linspace(-1e5,1e5,100)
    ax4.hist(ft1.fft3.imag.flatten(),histtype='step', bins=bins,label='imag')
    ax4.hist(ft1.fft3.real.flatten(),histtype='step', bins=bins,label='real')
    ax5.hist(ft2.fft3.imag.flatten(),histtype='step', bins=bins,label='imag')
    ax5.hist(ft2.fft3.real.flatten(),histtype='step', bins=bins,label='real')
    ax4.hist((np.abs(ft1.fft3)**2).flatten(), histtype='step',bins=bins,label='power')
    ax5.hist((np.abs(ft2.fft3)**2).flatten(), histtype='step',bins=bins,label='power')
    #ax4.hist(128**3*ft1.power.flatten(), histtype='step',bins=bins,label='power')
    #qqq=(np.abs(ft1.fft3)**2).flatten()/(128**3*ft1.power.flatten())
    #print(qqq)
    #ax5.hist(128**3*ft2.power.flatten(), histtype='step',bins=bins,label='power')
    ax5.legend(loc=0)
    ax4.legend(loc=0)
    ax4.set(yscale='log', xlim=[-5e5,5e5])
    ax5.set(yscale='log', xlim=[-5e5,5e5])

    ax6.hist(np.angle(ft1.fft3.flatten()))
    ax7.hist(np.angle(ft2.fft3.flatten()))



    fig.colorbar(p1,ax=ax3)
    fig.tight_layout()
    fig.savefig('plots_to_sort/cone_phases')


if 1:
    #Plot: projection, power vs radius, projection of FFT real
    fig,axes=plt.subplots(5,2,figsize=(12,12))
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ax4=axes[2][0];ax5=axes[2][1]
    ax6=axes[3][0];ax7=axes[3][1]
    ax8=axes[4][0];ax9=axes[4][1]

    p=ax0.imshow(np.log10(ft1.rho.sum(axis=0)))
    fig.colorbar(p,ax=ax0)
    p1=np.log10(ft2.rho.sum(axis=0))
    p1=np.roll(p1,p1.shape[0]//2,axis=0)
    p1=np.roll(p1,p1.shape[1]//2,axis=1)
    p=ax1.imshow(p1)
    fig.colorbar(p,ax=ax1)
    ok = rrr>0
    q1=np.log10(ft1.power[ok].flatten()).real
    q2=np.log10(ft2.power[ok].flatten()).real
    rr = np.log10(rrr[ok]).flatten()
    ext=dt.extents()
    ext(q1)
    ext(q2)
    binsy = np.linspace(ext.minmax[0],ext.minmax[1],64)
    binsx = np.linspace(rr.min(),rr.max(),64)
    bins=[binsx,binsy]
    pch.simple_phase(rr , q1, ax=ax2, bins=bins)
    ok = (ft2.power > 1e-15)*(rrr>0)
    pch.simple_phase( rr, q2, ax=ax3, bins=bins)

    
    binsx = np.linspace(theta.min(),theta.max(),64)
    pch.simple_phase(theta[ok].flatten(), q1, ax=ax6, bins=[binsx,binsy])
    pch.simple_phase(theta[ok].flatten(), q2, ax=ax7, bins=[binsx,binsy])

    binsx = np.linspace(phi.min(),theta.max(),64)
    pch.simple_phase(phi[ok].flatten(), q1, ax=ax8, bins=[binsx,binsy])
    pch.simple_phase(phi[ok].flatten(), q2, ax=ax9, bins=[binsx,binsy])


    q3= np.abs(ft1.power.real).sum(axis=0)
    q3=np.roll(q3,q3.shape[0]//2,axis=0)
    q3=np.roll(q3,q3.shape[0]//2,axis=1)
    q3=np.log(q3)
    ax4.imshow(q3)

    q3= np.abs(ft2.power.real).sum(axis=0)
    q3=np.roll(q3,q3.shape[0]//2,axis=0)
    q3=np.roll(q3,q3.shape[0]//2,axis=1)
    q3=np.log(q3)
    ax5.imshow(q3)


    fig.tight_layout()
    fig.savefig('plots_to_sort/why')

