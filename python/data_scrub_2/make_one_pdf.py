
from GL import *
from queb3 import powerlaw_fit as plfit

import simulation

import bucket

def make_pdf(simlist):
    
    fig,ax=plt.subplots(1,1)
    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]


        frame = this_sim.ann_frames[-1]

        pdf, bins,ds = bucket.things.get(sim,(None,None,None))
        if pdf is None:
            print('make pdf', sim)

            ds = this_sim.load_ds(frame)
            ad = ds.all_data()
            density = ad['density'].v.flatten()
            bins = np.geomspace(density.min(),density.max(),128)
            print('pdf')
            pdf,bins = np.histogram(density,bins=bins, density=True)
            bucket.things[sim]=pdf,bins,ds

        bc = 0.5*(bins[1:]+bins[:-1])

        ax.plot(bc, pdf,c=this_sim.color,linestyle=this_sim.linestyle)
    ax.set(xscale='log', yscale='log')
    outname="%s/pdf"%dl.plotdir
    fig.savefig(outname)
    print(outname)


