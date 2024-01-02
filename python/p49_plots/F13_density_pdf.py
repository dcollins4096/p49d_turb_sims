

from GL import *
from queb3 import powerlaw_fit as plfit

import simulation

import bucket

def plot_pdf(simlist):
    
    fig,ax=plt.subplots(1,1)
    for nsim,sim in enumerate(simlist):
        print('plot',sim)
        this_sim=simulation.corral[sim]

        this_sim.read_pdfs(['density'],pdf_prefix='pdf')

        bc = this_sim.cbins['density']
        pdf = this_sim.avg_pdf['density']
        ax.plot(bc, pdf,c=this_sim.color,linestyle=this_sim.linestyle)
    ax.set(xscale='log', yscale='log')
    outname="%s/pdf"%dl.plotdir
    fig.savefig(outname)
    print(outname)


