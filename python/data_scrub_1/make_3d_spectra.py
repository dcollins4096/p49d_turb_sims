print('import stuff')
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

print('word')

clobber=False
simlist=sim_colors.simlist
if len(sys.argv) > 1:
    simlist = [sys.argv[-1]]
else:
    simlist=sim_colors.simlist
#simlist = ['5_1']
for i,sim in enumerate(simlist):
    frames=sim_colors.framedict[sim]
    prefix='PREFIX'
    #sim_dir='/data/cb1/Projects/P49_EE_BB/%s'%sim
    #product_dir='/data/cb1/Projects/P49_EE_BB/Products/%s'%sim
    sim_dir='/scratch/00369/tg456484/Paper49/RUNNING/%s'%sim
    product_dir='/scratch/00369/tg456484/Paper49/RUNNING/Products/%s'%sim


    pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix, product_directory=product_dir)
    for frame in frames:
        pack.make_spectra(frame)
