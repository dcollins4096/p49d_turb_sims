from GL import *

import simulation as simulation
reload(simulation)
import compute_avg_quantities as comp_avg
reload(comp_avg)

def comp_all(simlist):
    for sim_name in simlist:
        this_sim = simulation.corral[sim_name]
        print(this_sim.all_frames)
        for frame in this_sim.ann_frames:
            if this_sim.code == 'Enzo':
                comp_avg.make_quan(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False, do_magnetic=this_sim.do_magnetic )
            elif this_sim.code == 'Athena':
                comp_avg.make_quan_athena(this_sim.data_location,frame,out_directory=this_sim.product_location,sim=sim_name, clobber=False, do_magnetic=this_sim.do_magnetic )


