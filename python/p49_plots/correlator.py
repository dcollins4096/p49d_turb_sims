
from GL import *

import simulation
reload(simulation)
import simulation_info.all_sims as all_sims

import do_correlation as dc
reload(dc)
#sim_list = all_sims.lists['suite1']
sim_list = ['6_half']

#eb.bin(sim_list)
dc.correlator(sim_list)


