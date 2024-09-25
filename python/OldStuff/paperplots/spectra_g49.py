
from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(sim_colors)
reload(queb3)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)


#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
#reload(read_stuff)  
spectra_dict = rs.spectra_dict
quan3 = rs.quan3

plotdir="%s/PigPen"%os.environ['HOME']
def smallest_string(val):
    if val > 0.6:
        return "%d"%int(val)
    else:
        return "%0.1f"%val
if 1:
    plt.close('all')
    #fig,ax = plt.subplots(1,3, figsize=(12,4))#, sharex=True,sharey=True,figsize=(12,4))
    #fig.subplots_adjust(wspace=0, hspace=0)
    #axlist=ax.flatten()
    fig,ax = plt.subplots(1,1)
    axlist=[ax]

    #fields=['avg_cltt','avg_clee','avg_clbb']
    fields=['avg_clee']
    for nf,field in enumerate(fields):
        for sim in sim_colors.simlist:
            ms,ma=sim_colors.vals_from_sim(sim)
            ms_string=smallest_string(ms)
            ma_string=smallest_string(ma)
            print(ms_string, ma_string)


            proj=rs.proj_dict['x'][sim]
            #label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            label="(%s, %s) %0.1f"%(ms_string,ma_string, spectra_dict['y'][sim].slopes[field])
            axlist[nf].plot(proj.lcent,   spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            axlist[nf].set_title(r'$C_\ell^{%s}$'%( field[-2:].upper()))

            #proj=rs.proj_dict['y'][sim]
            #label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            #axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            #axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    ax.plot( proj.fitrange, 100*(proj.fitrange/0.1)**(-2.5), c='k', label='Planck, -2.5')
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    print('poot')
    for a in axlist[:3]:
        a.legend(loc=0)
    for a in axlist:
        a.set_xlabel(r'$\ell/\ell_{max}$')
    for a in [axlist[0]]:#,axlist[3]]:
        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/alpha_TEB.pdf'%plotdir)

