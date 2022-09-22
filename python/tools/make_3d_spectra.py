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

clobber=False
simlist=sim_colors.simlist
if len(sys.argv) > 1:
    simlist = [sys.argv[-1]]
else:
    simlist=['5_2', '5_3']#,'5_half','5_1']

for i,sim in enumerate(simlist):
    frames=sim_colors.framelist[i]
    prefix='PREFIX'
    sim_dir='/data/cb1/Projects/P49_EE_BB/%s'%sim

    pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix)
    for frame in frames:
        pack.make_spectra(frame)
