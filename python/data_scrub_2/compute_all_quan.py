from GL import *

import simulation as sim
reload(sim)
import compute_avg_quantities as comp_avg
reload(comp_avg)

def comp_all(sim_name):
    this_sim = sim.corral[sim_name]
    for frame in this_sim.framelist:
        comp_avg.make_quan(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name )

