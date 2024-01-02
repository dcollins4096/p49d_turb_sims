


from GL import *

import simulation
reload(simulation)
reload(sim_colors)
def plot_summary(simlist, LOS = 'x',ax=None):
    
    savefig=False
    if ax is None:
        fig,ax=plt.subplots(1,1)
        savefig=True
    plotdir =  "/home/dccollins/PigPen"

    for sim in simlist:
        this_sim = simulation.corral[sim]
        this_sim.load()
        simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)
        #the_x=spectra_dict[LOS][sim].amps['avg_clbb']/spectra_dict[LOS][sim].amps['avg_clee']
        #the_y=spectra_dict[LOS][sim].slopes['avg_clee']
        the_x = this_sim.ampsA['ClBB'+LOS]/this_sim.ampsA['ClEE'+LOS]
        the_y = this_sim.slopesA['ClEE'+LOS]
        kwargs = {"c":[this_sim.color], "marker":this_sim.marker, "s":this_sim.marker_size*100}
        plot=ax.scatter( the_x, the_y, **kwargs)
        ax.set(xlabel=r'$A_{BB}/A_{EE}$', ylabel=r'$\alpha_{EE}$')

    colorbar=plt.colorbar(sim_colors.cbar,ax=ax)
    colorbar.set_label(sim_colors.mach_label)
    ax.axhline(sim_colors.planck_E_slope,c= [0.5]*4)
    ax.axvline(sim_colors.planck_ratio,c= [0.5]*4)
    if savefig:
        fig.tight_layout()
        outname = '%s/summary.pdf'%dl.plotdir
        fig.savefig(outname)
        print(outname)
