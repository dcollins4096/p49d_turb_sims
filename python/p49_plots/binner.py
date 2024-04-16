
from GL import *

import simulation
reload(simulation)
import simulation_info.all_sims as all_sims

import equal_bins as eb
reload(eb)
#sim_list = all_sims.lists['suite1']
sim_list = ['6_half']

#eb.bin(sim_list)
eb.parabin(sim_list)


