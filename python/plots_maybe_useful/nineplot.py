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

plot_dir =  "/home/dccollins/PigPen"
plotdir =  "/home/dccollins/PigPen"

#
# Primitive Slope vs TEB slope
#
LOS = 'x'
if 1:
    plt.close('all')
    fig,ax = plt.subplots(3,3,figsize=(12,4))#, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nfx,fieldx in enumerate(['avg_d','avg_v','avg_h']):
        for nfy,fieldy in enumerate(['avg_cltt','avg_clee','avg_clbb']):
            for sim in sim_colors.simlist:
                #label="%s %0.2f"%(sim, spectra_dict['x'][sim].slopes[field])
                #print(label)
                #kwargs={'c':sim_colors.color[sim], 'linestyle':sim_colors.linestyle[sim], 
                #        'label':label}
                ax[nfy][nfx].scatter( spectra_dict[LOS][sim].slopes[fieldx],
                                    spectra_dict[LOS][sim].slopes[fieldy],
                                   c=sim_colors.color[sim], marker=sim_colors.marker[sim])
                #kwargs.pop('label')
                #pline(spectra_dict['x'][sim], axlist[nf], field, c='k')

            
    for a in axlist:
        dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
    ax[2][0].set_xlabel(r'$\alpha_\rho$')
    ax[2][1].set_xlabel(r'$\alpha_v$')
    ax[2][2].set_xlabel(r'$\alpha_H$')
    ax[0][0].set_ylabel(r'$\alpha_{\rm{TT}}$')
    ax[1][0].set_ylabel(r'$\alpha_{\rm{EE}}$')
    ax[2][0].set_ylabel(r'$\alpha_{\rm{BB}}$')
    ax[0][0].set_title('take1')
    ax[0][-1].yaxis.tick_right()
    ax[1][-1].yaxis.tick_right()
    ax[2][-1].yaxis.tick_right()
    ax[0][1].yaxis.set_ticklabels([])
    ax[1][1].yaxis.set_ticklabels([])
    ax[2][1].yaxis.set_ticklabels([])
    for aaa in ax.flatten():
        aaa.grid()


    outname = '%s/slope_TTEEBB_vs_rho_v_H.pdf'%plotdir
    fig.savefig(outname)
    print(outname)
