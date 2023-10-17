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
    for sim in simlist:
        print("All spec",sim)
        this_sim = simulation.corral[sim]
        #sim_dir = "/scratch/00369/tg456484/Paper49/%s"%sim
        #product_dir = "/scratch/00369/tg456484/Paper49/Products/%s"%sim

#a thing that describes the simulation
        prefix = this_sim.name
        pack = queb3.simulation_package( directory=this_sim.data_location,frames=this_sim.all_frames,prefix=prefix, 
                                        product_directory=this_sim.product_location, simname=sim)
#produce all QUEB products.
        #pack.EBall()
        pack.ebspec()
