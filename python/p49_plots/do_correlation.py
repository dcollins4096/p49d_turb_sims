
from GL import *

import simulation
reload(simulation)

def correlator(sim_list):


    for sim in sim_list:
        this_sim = simulation.corral[sim]
        this_sim.load()

        #ann_frames is the list of analysis frames for this simulation.
        frames = this_sim.ann_frames[-1:]

        b = this_sim.get_field('b', frames[0], ax='y')
        t = this_sim.get_field('d', frames[0], ax='y')
        e = this_sim.get_field('d', frames[0], ax='y')

