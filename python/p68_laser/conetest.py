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
L=-0.5; R=0.5; Dx=1/N
kI = np.fft.fftfreq(N)*N
kx,ky,kz=np.meshgrid(kI,kI,kI)
rrr = np.sqrt(kx**2+ky**2+kz**2)
Ahat = np.zeros_like(rrr*1j)
ok = rrr>0
Ahat[ok]=rrr[ok]**-1.5
if 0:
    phi = np.random.random(Ahat.size)
    phi.shape = Ahat.shape
    Ahatmag = np.abs(Ahat)
    Ahat = Ahatmag*np.cos(phi)+Ahatmag*np.sin(phi)*1j
H = Ahat.shape[0]//2
Ahat = Ahat[:,:,:H+1]
Q = np.fft.irfftn(Ahat)
#Aback = np.fft.ifftn(Q).real
#fig,axes=plt.subplots(2,2)
#ax0=axes[0][0]

#ax0.imshow(Q.sum(axis=0))
#fig.savefig('plots_to_sort/herp')

#
# kludge
#

def fitter(X,Y):
    ok=(X>0)*(Y>0)
    pfit = np.polyfit(np.log10(X[ok]),np.log10(Y[ok]),1)
    return pfit

if 1:
    shell = bt.shell_average(Ahat.real)
    TheX_x = shell.kd[1:]*128
    TheY_x = shell.power_1d[1:]/shell.Nzones[1:]
    ok_x = (TheX_x>0)*(TheY_x>0)
    pfit_x = fitter(TheX_x,TheY_x)
    print('I get',pfit_x)

if 1:
    #FFT = np.fft.fftn(Q)
    #power = FFT*np.conj(FFT)
    #power = (Ahat*Ahat.conj()).real
    power = Ahat.real
    kbins = np.sqrt( np.unique(kI**2))
    kbins.sort()
    kcen = 0.5*(kbins[1:]+kbins[:-1])

    kvals = rrr#[:,:,:H+1]
    TheY, bin_edge, count = scipy.stats.binned_statistic(kvals.flatten(), power.flatten(), bins=kbins, statistic='mean')

    TheX = kcen
    ok = (TheX>0)*(TheY>0)
    LogX=np.log10(TheX[ok])
    pfit = fitter(TheX,TheY)
    print('I get',pfit)
    fig,axes=plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]
    ax0.scatter( kvals.flatten(), power.flatten(), s=0.1)
    ext=dt.extents()
    ext(power,positive=True)
    extX=dt.extents()
    extX(kvals,positive=True)
    TheY2=10**(pfit[0]*LogX+pfit[1])
    #ax0.plot( 10**LogX, TheY2,c='g')
    #ax0.plot(TheX, TheY,c='g')
    ext(TheY2)
    ax0.set(yscale='log',xscale='log')#,ylim=ext.minmax, xlim=extX.minmax)
    ax0.plot( shell.kd*128, shell.power_1d/shell.Nzones, c='r')
    ax0.plot(TheX_x[ok_x],10**(pfit_x[0]*np.log10(TheX_x[ok_x])+pfit_x[1]),c='g')


    ax1.plot(TheX_x,np.sqrt(shell.Nzones[1:]))

    fig.savefig('plots_to_sort/law')

#if 'ftool' not in dir():
if 0:
    ftool=bt.fft_tool(Q)
    ftool.do3()
    ftool.do2(projax=0)
if 0:
    bt.plot_fft(ftool, outname ='plots_to_sort/cone')
if 0:
    slabimg.plot_fft(ftool, outname = "plots_to_sort/cone_fft")
    if 'slab' not in dir():
        slab=slabimg.slab(ftool, projax=0)
    slabimg.number_checker(ftool, outname = "plots_to_sort/cone_number")
