from GL import *

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
#the import does the reading.
#import read_stuff as rs
#reload(read_stuff)  
plotdir = dl.plotdir

import simulation
import simulation_info.all_sims

def nom(simlist):
    plt.close('all')
    fig,ax = plt.subplots(1,1,sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for i,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]
        this_sim.load()
        #ax.scatter( raq.quan3[sim]['maavg'], raq.quan3[sim]['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim],s=60)
        ax.scatter( this_sim.Ma_mean, this_sim.Ms_mean,c=this_sim.color, marker=this_sim.marker,s=60)
    dt.axbonk(ax,xlabel=r'$M_{\rm{A}}$',ylabel=r'$M_{\rm{S}}$')#,xlim=[0.45,2.1],ylim=[0.45,3.2])
    fig.savefig('%s/point_legend_measured.pdf'%plotdir)
            

