

from starter1 import *
import yt
from downsample import volavg
import fourier_tools_py3.fourier_filter as Filter
import tools.davetools as dt

def plot_fft(self, outname="./ploot"):
    ext=dt.extents()
    ext(self.power.real[self.power.real>0])

    if 0:
        shape = self.power.real.shape
        xx,yy,zz = np.mgrid[0:shape[0]:1, 0:shape[1]:1, 0:shape[2]:1]+1
        ext=dt.extents()
        ext(xx)
    for S in [0]:
        print("Slab",S)
        fig,axes=plt.subplots(2,2)
        axlist=axes.flatten()

        fx={}
        for d in range(3):
            zero = [slice(None),slice(None),slice(None)]
            zero[d]=S

            #thing=self.power.sum(axis=d).real
            thing = self.power[tuple(zero)].real
            thing=np.roll(thing, thing.shape[0]//2, axis=0)
            thing=np.roll(thing, thing.shape[1]//2, axis=1)
            if 0:
                thing = yy[tuple(zero)]
                print(thing.shape)
                print(thing)
            fx[d]= thing
            
            #ext(fx[d][ fx[d]>0] )



        norm = mpl.colors.LogNorm(vmin=ext.minmax[0],vmax=ext.minmax[1])
        #norm = mpl.colors.Normalize(vmin=ext.minmax[0],vmax=ext.minmax[1])
        print(norm)


        axlist[0].imshow(fx[0],norm=norm, origin='lower', interpolation='nearest')
        axlist[1].imshow(fx[1],norm=norm, origin='lower', interpolation='nearest')
        axlist[2].imshow(fx[2],norm=norm, origin='lower', interpolation='nearest')
        axlist[0].set(xlabel='Z', ylabel='Y')
        axlist[1].set(xlabel='Z', ylabel='X')
        axlist[2].set(xlabel='Y', ylabel='X')

        fig.tight_layout()
        fig.savefig(outname + "%04d"%S)
        print("save")

def slab(self,projax=0):
    print('slab')
    self.fft3 = np.fft.fftn( self.rho )
    self.power=self.fft3*np.conjugate(self.fft3)
    self.power/=self.power.size
    ff = Filter.FourierFilter(self.power)
    nx,ny,nz=self.power.shape
    print(nx,ny,nz)
    pdb.set_trace()
    self.slab_spectrum = []
    my_coord = ff._kk[projax]


    self.slff=ff
    from scipy import stats
    for the_k in my_coord.flatten():
        print(" do",the_k)
        filtered_by_k = (my_coord==the_k)*self.power
        #spectrum, bin_edge, bin_num = stats.binned_statistic( self.slff.kradius, filtered_by_k, statistic='sum',bins=ff.bins)
        #self.slab_spectrum.append(spectrum)
        self.slab_spectrum.append(  np.array([filtered_by_k[ff.get_shell(bin)].sum() for bin in range(ff.nx)]))

    #self.power_1d3 /= self.rho.size
def number_checker(ftool,outname=None,axin=None,label=None):
    sigma_2 = (ftool.rho**2).sum()
    sigma_k = (ftool.power_1d3.real).sum()
    sigP_2 = (ftool.rho2**2).sum()
    sigP_k = (ftool.power_1d2.real).sum()
    sigPPk = (ftool.power_1d2.real*ftool.k2d).sum()

    slab = np.array(ftool.slab_spectrum)

    #sigPPk = sigma_k #just to check.
    if 1:
        #things that need to work
        print("sigma r",sigma_2)
        print("sigma k",sigma_k)
        print("sigmP r",sigP_2)
        print("sigmP k",sigP_k)
    if 1:
        sig2_guess = sigP_2*sigma_k/sigP_k
        print("should work", sig2_guess/sigma_2)
        sig2_guess = sigP_2*sigPPk/sigP_k
        print("brunt ratio sigma^2 brunt/sigma^2 actual", sig2_guess/sigma_2)
        print("One?",slab.sum()/sigma_k)
        if axin is None:
            fig,ax=plt.subplots(1,1)
        else:
            ax=axin
        q=slab.sum(axis=0)
        #q=np.cumsum(q)
        ax.plot(ftool.k2d,q/q[20],label=label)
        ax.set(yscale='log')
        ax.set_xscale('symlog')
        if axin is None:
            fig.savefig(outname)


def plot(ftool,outname="./plot"):


    fig, axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
        
    sl=slice(1,None)
    q=ftool.power_1d2.real[sl]
    ax0.plot(ftool.k2d[sl], q/q[20],c='r', label='P2d')
    q=(ftool.slab_spectrum[0].real*128*10)[sl]
    ax0.plot(ftool.k2d[sl], q/q[20], label='P3d(k=0)')
    q=(ftool.power_1d3.real/ftool.k3d)[sl]
    ax0.plot(ftool.k3d[sl], q/q[20], c='g', label='p3d/k')
    q=ftool.power_1d3.real[sl]
    ax0.plot(ftool.k3d[sl], q/q[20], 'g--', label='p3d')
    ax0.legend(loc=0)
    slab = np.array(ftool.slab_spectrum)
    ax1.plot(ftool.k3d[sl], ftool.power_1d3.real[sl], c='g',label='p3d')
    ax1.plot(ftool.k3d[sl], slab.sum(axis=0).real[sl], label=r'$\sum P3d(k=k)$')
    ax1.legend(loc=0)
    ax2.plot(ftool.k3d[sl],ftool.power_1d3.real[sl], c='g',label='p3d')
    ax2.plot(ftool.k3d[sl],slab[1:,:].sum(axis=0).real[sl], c='r',label='P3d(k>0)')
    ax2.plot(ftool.k3d[sl],slab[0:1,:].sum(axis=0).real[sl], label='P3d(k==0)')
    ax2.legend(loc=0)
    ax3.plot(ftool.k3d[sl], ftool.power_1d3.real[sl], c='g',label='p3d')
    for kz in range(slab.shape[0]):
        ok = slab[kz,:]>0
        ok[0]=False
        c=[0.5]*4
        if kz == 0:
            c = 'k'
        ax3.plot(ftool.k3d[ok], slab[kz,ok].real,c=c)
    ax0.set(xscale='log',yscale='log')
    ax1.set(xscale='log',yscale='log')
    ax2.set(xscale='log',yscale='log')
    ax3.set(xscale='log',yscale='log')


    fig.savefig(outname)
    print("wrote %s"%outname)

    plt.close(fig)
