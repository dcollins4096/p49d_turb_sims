"""
Make many queb data products.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
import queb3
reload(queb3)
import simulation

def make_all(simlist):
    print('simlist',simlist)
    for sim in simlist:
        print(sim)
        this_sim = simulation.corral[sim]
        #sim_dir = "/scratch/00369/tg456484/Paper49/%s"%sim
        #product_dir = "/scratch/00369/tg456484/Paper49/Products/%s"%sim
        if this_sim.B_nom > 0:
            do_magnetic=True
        else:
            do_magnetic=False

#a thing that describes the simulation

        for frame in this_sim.ann_frames:
            prefix = this_sim.name
            pack = queb3.simulation_package( directory=this_sim.data_location,frames=[frame],prefix=prefix, 
                                            product_directory=this_sim.product_location, simname=sim, code=this_sim.code)
#produce all QUEB products.
            pack.EBall(do_magnetic=do_magnetic)
            del pack
            import gc
            gc.collect()
