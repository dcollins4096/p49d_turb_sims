from GL import *

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
#the import does the reading.
#import read_stuff as rs
#reload(read_stuff)  
plotdir = dl.plotdir

import simulation
reload(simulation)
import simulation_info.all_sims
sim_colors.cmap='gist_rainbow'
reload(sim_colors)

def nom(simlist):
    plt.close('all')
    fig,ax = plt.subplots(1,1,sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for i,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]
        this_sim.load()
        simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)
        #ax.scatter( raq.quan3[sim]['maavg'], raq.quan3[sim]['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim],s=60)
        ax.scatter( this_sim.Ma_mean, this_sim.Ms_mean,c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
    colorbar=fig.colorbar(sim_colors.cbar,ax=ax)
    colorbar.set_label(sim_colors.mach_label)
    #pdb.set_trace()
    print('args')
    
    dt.axbonk(ax,xlabel=sim_colors.alf_mach_label,ylabel=sim_colors.mach_label)
    outname='%s/point_legend_measured.pdf'%plotdir
    fig.savefig(outname)
    print(outname)
            

