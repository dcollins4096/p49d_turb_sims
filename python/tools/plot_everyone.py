


from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(queb3)
reload(sim_colors)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)
plot_dir="/home/dccollins/PigPen"
clobber=False
simlist=sim_colors.simlist
#simlist=['5_half','5_1','5_2', '5_3']
shortprefix="time"

import read_avg_spectra as ras

for axes in ['x','y']:
    for i,sim in enumerate(simlist):
        print('SIM',sim)
        frames=sim_colors.framelist[i]
        sim_dir = "/data/cb1/Projects/P49_EE_BB/%s"%sim
        product_dir = "/data/cb1/Projects/P49_EE_BB/Products/%s"%sim
        pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix='test', product_directory=product_dir)
        fig,ax=plt.subplots(1,2)

        for frame in frames:
            ax[0].clear()
            ax[1].clear()
            proj=pack.read_queb(frame=frame,ax=axes,bin_style='dx1')
            proj.compute_harmonic_products()
            norm = mpl.colors.LogNorm(vmin=0.5,vmax=2)
            ax[0].imshow(proj['T'], norm=norm)
            spectra =  ras.spectra_dict[axes][sim]


            ax[1].plot(spectra.lcent, spectra.spectra['avg_cltt'], 'k')

            ax[1].plot(spectra.lcent, proj['ClTT'])
            ax[1].set(xscale='log',yscale='log', title='ClTT %s %04d'%(sim,frame))
            outname='%s/eyes_%s_%s_%04d'%(plot_dir,sim,axes,frame)
            fig.savefig(outname)
            print(outname)

