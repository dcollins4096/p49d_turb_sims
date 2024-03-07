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

ax0.plot(TNA['r120_t1'].sum(axis=0),c='r')
ax0.plot(TNA['r120_t2'].sum(axis=0),c='g')
ax0.set_title('r120')

ax1.plot(TNA['r60_t1'].sum(axis=0),c='r')
ax1.plot(TNA['r60_t2'].sum(axis=0),c='g')
ax1.set_title('r60')
ax2.plot(TNA['r0_t1'].sum(axis=0),c='r')
ax2.plot(TNA['r0_t2'].sum(axis=0),c='g')
ax2.set_title('r0')
fig.savefig('plots_to_sort/shock_fronts')
