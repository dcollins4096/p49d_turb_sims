


from GL import *

import simulation

def plot_summary(simlist, LOS = 'x'):
    
    fig,ax=plt.subplots(1,1)
    plotdir =  "/home/dccollins/PigPen"

    for sim in simlist:
        this_sim = simulation.corral[sim]
        this_sim.load()
        #the_x=spectra_dict[LOS][sim].amps['avg_clbb']/spectra_dict[LOS][sim].amps['avg_clee']
        #the_y=spectra_dict[LOS][sim].slopes['avg_clee']
        the_x = this_sim.ampsA['ClBB'+LOS]/this_sim.ampsA['ClEE'+LOS]
        the_y = this_sim.slopesA['ClEE'+LOS]
        kwargs = {"c":this_sim.color, "marker":this_sim.marker}
        ax.scatter( the_x, the_y, **kwargs)
        ax.set(xlabel=r'$A^{BB}/A^{EE}$', ylabel=r'$\alpha_{EE}$')

    ax.axhline(-2.45,c= [0.5]*4)
    ax.axvline(0.5,c= [0.5]*4)
    fig.tight_layout()
    outname = '%s/summary.pdf'%dl.plotdir
    fig.savefig(outname)
    print(outname)
