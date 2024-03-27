from GL import *

import simulation
reload(simulation)
import simulation_info.all_sims as all_sims


sim_list = all_sims.lists['suite1']


if 1:
    import F15_driving as F15
    reload(F15)
    sim_list=['half_half','1_2','6_half','6_2']
#F15.plot_accel(sim_list)
    F15.plot_driving(sim_list)

