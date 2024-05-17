
from dtools.starter1 import *
import dtools.davetools as dt
import power_spectrum as ps
reload(ps)

class noisy():
    def __init__(self, arr):
        self.arr=arr
        self.ps = ps.powerspectrum(self.arr)
    def kpicker(self,kmin=None,kmax=None, nmin=None,nmax=None):
        k1 = self.ps.kcen.min()
        k2 = self.ps.kcen.max()
        if nmin is not None:
            k1 = self.ps.kcen[nmin]
        if nmax is not None:
            k2 = self.ps.kcen[nmax]
        if kmin is not None:
            k1 = kmin
        if kmax is not None:
            k2 = kmax
        return k1, k2
    def denoise_1(self,kmin=None,kmax=None, nmin=None,nmax=None):
            #rho=preshock_region[shot]
            Nhat=self.ps.Nhat+0
            rhohat = np.abs(Nhat)**2
            k1,k2=self.kpicker(kmin,kmax,nmin,nmax)
            print('K1 K2',k1,k2)
            Nhat[ self.ps.k < k1]=0
            Nhat[ self.ps.k > k2]=0
            rhoback = np.fft.ifftn(Nhat)
            self.Nhatfiltered=Nhat
            self.rhoback=rhoback
            return self.rhoback

    def plot_denoise(self,kmin=None,kmax=None,nmin=None,nmax=None,outname='plot'):
        self.denoised = self.denoise_1(kmin=kmin,kmax=kmax,nmin=nmin,nmax=nmax)
        fig,axes=plt.subplots(2,2, figsize=(12,12))
        ax0=axes[0][0];ax1=axes[0][1]#;ax2=axes[0][2]
        ax3=axes[1][0];ax4=axes[1][1]#;ax5=axes[1][2]
        #ax0.hist( self.arr.flatten(),histtype='step')
        #ax0.set(xlabel='Ival',ylabel='N')

        norm = mpl.colors.LogNorm(vmin=self.ps.rhohat[self.ps.rhohat>0].min(), vmax=self.ps.rhohat.max())
        ax0.imshow(self.ps.rhohat,norm=norm)

        #ax3.plot(self.ps.kcen,self.ps.power, marker='*')
        knorm=self.ps.kcen.min()
        ax3.plot(self.ps.kcen/knorm,self.ps.power, marker='*')
        ax3.set(xscale='log',yscale='log')

        k1,k2=self.kpicker(kmin,kmax,nmin,nmax)
        ax3.axvline(k1/knorm)
        ax3.axvline(k2/knorm)

        A1 = self.arr
        A2 = self.rhoback.real
        ext = dt.extents()
        ext(A1)
        ext(A2)
        norm = mpl.colors.Normalize(ext.minmax[0],ext.minmax[1])
        plot=ax0.imshow(A1,norm=norm)
        cb=fig.colorbar(plot,ax=ax0)
        #cb.set_title('IVAL')
        plot=ax1.imshow(A2,norm=norm)
        cb2=fig.colorbar(plot, ax=ax1)
        #ax5.imshow(A1-A2)
        #ax5.imshow(np.log(np.abs(self.Nhatfiltered)))
        ax4.plot(A1[:,250])
        ax4.plot(A2[:,250])











        print(outname)
        fig.tight_layout()
        fig.savefig(outname)




