from dtools.starter1 import *


if 'TNA' not in dir():
    fname = 'p68_laser/TRIM_ALIGN.h5'
    fptr=h5py.File(fname,'r')
    TNA = {}
    for field in fptr:
        TNA[field]=fptr[field][()]
    fptr.close()

fig,axes=plt.subplots(3,1)
ax0=axes[0];ax1=axes[1];ax2=axes[2]

def doq(name, ax, **kwargs):
    Q = TNA[name]
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
        x = np.arange(Qbar.size)
        #ax.plot(Qbar, **kwargs)
        #ax.plot(Qbar+std, linewidth=0.2, **kwargs)
        #ax.plot(Qbar-std, linewidth=0.2, **kwargs)
        ax.plot(std,**kwargs)

doq('r120_t1',ax0,c='r')
doq('r120_t2',ax0,c='g')
ax0.set_title('r120')

doq('r60_t1',ax1,c='r')
doq('r60_t2',ax1,c='g')
ax1.set_title('r60')
doq('r0_t1',ax2,c='r')
doq('r0_t2',ax2,c='g')
ax2.set_title('r0')
fig.savefig('plots_to_sort/shock_fronts')
