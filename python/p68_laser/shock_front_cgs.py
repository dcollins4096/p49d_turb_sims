from dtools.starter1 import *

import regions
reload(regions)
import physical_values as phys
reload(phys)
if 'TNA' not in dir():
    fname = 'p68_laser/TRIM_ALIGN.h5'
    fptr=h5py.File(fname,'r')
    TNA = {}
    for field in fptr:
        TNA[field]=fptr[field][()]
    fptr.close()


def doq(name, ax,ax2=None, **kwargs):
    image = TNA[name]
    ps = regions.preshock_region[name]
    I0 = phys.compute_I0(ps).mean()
    zero_region = regions.zero_region[name]
    zero_value = 2*zero_region.min()
    image_sort = copy.copy(image.flatten())
    image_sort.sort()
    delta = image_sort[1]-image_sort[0] 
    zero_value = image_sort[0]-delta
    Q = phys.compute_rho(image, I0, zero=zero_value)
    print("Negative values:", (Q<=0).sum())

    if 0:
        for nx,ix in enumerate(range(0,Q.shape[0],10)):
            print(name,ix)
            ax.plot(Q[ix,:],linewidth=0.1,**kwargs)
        ax.plot( Q.mean(axis=0),**kwargs)
    if 0:
        std = Q.std(axis=0)
        Qbar = Q.mean(axis=0)
        x = np.arange(Qbar.size)
        ax.errorbar(x,Qbar, yerr=std)
    if 1:
        sl = slice(200,300)
        std = Q[sl,:].std(axis=0)
        Qbar = Q[sl,:].mean(axis=0)
        ax.plot(Q[sl,:].transpose(), linewidth=0.2,c=[0.5]*4)
        ax.plot(Qbar, **kwargs)
        #ax.plot(Q[600,:],c='k')
        #ax.plot(Qbar+std, linewidth=0.2, **kwargs)
        #x = np.arange(Qbar.size)
        #ax.plot(Qbar-std, linewidth=0.2, **kwargs)
        #ax.plot(std,**kwargs)
        ax.axhline(0)
        #ax.set(yscale='log')
    if ax2:
        xax = np.arange(0,Q.shape[0])
        ax2.imshow(Q)
        for line in xax[sl]:
            ax2.axhline(line,c=[0.5,0.5,0.5,0.1])

fig,axes=plt.subplots(3,3)
ax0=axes[0][0];ax1=axes[0][1];ax2=axes[0][2]
ax3=axes[1][0];ax4=axes[1][1];ax5=axes[1][2]
ax6=axes[2][0];ax7=axes[2][1];ax8=axes[2][2]
doq('r0_t1',ax0,ax1,c='r')
doq('r0_t2',ax0,ax2,c='g')
doq('r60_t1',ax3,ax4,c='r')
doq('r60_t2',ax3,ax5,c='g')
doq('r120_t1',ax6,ax7,c='r')
doq('r120_t2',ax6,ax8,c='g')
fig.tight_layout()
fig.savefig('plots_to_sort/shock_fronts')

if 0:
    fig,axes_sq=plt.subplots(3,2)
    axes=axes_sq.flatten()
    ax0=axes[0];ax1=axes[1];ax2=axes[2]
    ax3=axes[3];ax4=axes[4];ax5=axes[5]
    doq('r120_t1',ax0,ax3,c='r')
    doq('r120_t2',ax0,c='g')
    ax0.set_title('r120')

    doq('r60_t1',ax1,c='r')
    doq('r60_t2',ax1,c='g')
    ax1.set_title('r60')
    doq('r0_t1',ax2,c='r')
    doq('r0_t2',ax2,c='g')
    ax2.set_title('r0')
    fig.savefig('plots_to_sort/shock_fronts')
