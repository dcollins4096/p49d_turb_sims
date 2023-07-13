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

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
#the import does the reading.
#import read_stuff as rs
#reload(read_stuff)  
plotdir = dl.plotdir

import simulation

if 1:
    plt.close('all')
    fig,ax = plt.subplots(1,1,sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for i,sim in enumerate(sim_colors.simlist):
        this_sim=simulation.corral[sim]
        this_sim.read_avg_quan()
        #ax.scatter( raq.quan3[sim]['maavg'], raq.quan3[sim]['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim],s=60)
        ax.scatter( this_sim.quan3['maavg'], this_sim.quan3['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim],s=60)
    dt.axbonk(ax,xlabel=r'$M_{\rm{A}}$',ylabel=r'$M_{\rm{S}}$')#,xlim=[0.45,2.1],ylim=[0.45,3.2])
    fig.savefig('%s/point_legend_measured.pdf'%plotdir)
            

