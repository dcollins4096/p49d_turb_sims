"""
Make many queb data products.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
import sim_colors
import queb3
reload(queb3)
sim = sys.argv[1]
frames=sim_colors.framedict[sim]
#sim_dir = "/scratch/00369/tg456484/Paper49/%s"%sim
#product_dir = "/scratch/00369/tg456484/Paper49/Products/%s"%sim
sim_dir = sim_colors.cloudbreak_base + "/" + sim
product_dir = sim_colors.cloudbreak_base + "/Products/" + sim
prefix=sim

#a thing that describes the simulation
pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix, product_directory=product_dir)
#produce all QUEB products.
pack.EBall()
