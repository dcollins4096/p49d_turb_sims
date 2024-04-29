
from dtools.starter1 import *
import dtools.davetools as dt
import power_spectrum as ps
import get_foam as gf
reload(gf)

import preshock_regions as preshock

shot='r60_t1'
arr =    preshock.preshock_region[shot]
arrhat = ps.powerspectrum(arr)

fig,ax=plt.subplots(2,2, figsize=(12,12))
axes=ax.flatten()

def b(qqq):
    a = np.zeros( nar(qqq.shape)//2+nar([1,1]))
    nx,ny=qqq.shape
    a[1:nx//2+1,1:ny//2+1]= qqq[:qqq.shape[0]//2,:qqq.shape[1]//2]
    return a

Nhatreal=b(np.abs(arrhat.Nhat))
Nhatimag=np.cos(b(np.angle(arrhat.Nhat)))
#Nhatreal=b(arrhat.Nhat.real)
#Nhatimag=b(arrhat.Nhat.imag)
print(Nhatreal.shape)
print(arrhat.Nhat.shape)
norm = mpl.colors.SymLogNorm(vmin=Nhatreal.min(),vmax=Nhatreal.max(),linthresh=1)
p=axes[0].imshow(Nhatreal,norm=norm, origin='lower')
fig.colorbar(p,ax=axes[0])
#norm = mpl.colors.SymLogNorm(vmin=Nhatimag.min(),vmax=Nhatimag.max(),linthresh=1)
norm = mpl.colors.Normalize(vmin=Nhatimag.min(),vmax=Nhatimag.max())
p=axes[1].imshow(Nhatimag,norm=norm, origin='lower', cmap='twilight')
fig.colorbar(p,ax=axes[1])
p=axes[2].imshow(arrhat.rho)
fig.colorbar(p,ax=axes[2])
xlim=[1,Nhatreal.shape[0]-1]
ylim=[1,Nhatreal.shape[1]-1]
axes[0].set(xscale='log',yscale='log', xlim=xlim,ylim=ylim)
axes[1].set(xscale='log',yscale='log', xlim=xlim,ylim=ylim)
axes[3].hist( Nhatimag.flatten(), histtype='step')
fig.savefig('plots_to_sort/fft_only_%s'%shot)

fig,axes=plt.subplots(3,3)
axlist=axes.flatten()

norm = mpl.colors.Normalize(vmin=arr.min(),vmax=arr.max())
p=axlist[0].imshow(arrhat.rho,norm=norm)
fig.colorbar(p,ax=axlist[0])
nhat=np.abs(arrhat.Nhat)
norm2 = mpl.colors.LogNorm(vmin=nhat.min(),vmax=nhat.max())
axlist[1].imshow(nhat,norm=norm2)

p=axlist[2].plot( arrhat.kcen, arrhat.power)
axlist[2].set(xscale='linear',yscale='log')
#axlist[2].set(xscale='log',yscale='log')

def band_pass(arrhat,kmin,kmax, ax,norm=None):
    dk = arrhat.k[0,1] #first point.  do better later.
    ok = (arrhat.k >= kmin*dk)*(arrhat.k<=kmax*dk)
    print('fraction %0.2e'%(ok.sum()/ok.size))
    filt_hat = np.zeros_like(arrhat.Nhat)
    filt_hat[ok] = arrhat.Nhat[ok]
    filt_real = np.fft.ifftn(filt_hat)
    if norm is None:
        norm = mpl.colors.Normalize(vmin=filt_real.real.min(),vmax=filt_real.real.max())
    p=ax.imshow(filt_real.real,norm=norm)
    print("Sum",filt_real.real.sum())
    fig.colorbar(p,ax=ax)

band_pass(arrhat, 0,0.9,axlist[3],norm=norm)
axlist[3].set(title='k=0')
band_pass(arrhat, 1,2-1e-9,axlist[4],norm=None)
axlist[4].set(title='k=1')
band_pass(arrhat, 2,3-1e-9,axlist[5],norm=None)
axlist[5].set(title='k=2')
band_pass(arrhat, 3,4-1e-9,axlist[6],norm=None)
axlist[6].set(title='k=3')
band_pass(arrhat, 4,5-1e-9,axlist[7],norm=None)
axlist[7].set(title='k=4')
band_pass(arrhat, 5,1e7,axlist[8],norm=None)
axlist[8].set(title='k>5')


fig.tight_layout()
fig.savefig('plots_to_sort/fft_games_%s'%shot)
