

from GL import *
from downsample import volavg

plot_dir = dl.plotdir

import simulation
reload(simulation)
import simulation_info.all_sims
import brunt_tools as bt
reload(bt)

def make(sim_list):

    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        frame = this_sim.ann_frames[-1]
        rho = this_sim.load_small_rho(frame)
