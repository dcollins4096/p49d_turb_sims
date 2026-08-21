
from GL import *
import simulation
reload(simulation)
import sim_colors
import re

list_and_frames = [["xi_0.25_mach112", 300],
                   ["xi_0.25_mach12", 300],
                   ["xi_0.25_mach16", 300],
                   ["xi_0.25_mach31.1", 300],
                   ["xi_0.25_mach5.1", 300],
                   ["xi_0.25_mach59", 300],
                   ["xi_0.25_mach6.8", 300],
                   ["xi_0.5_mach1", 301],
                   ["xi_0.5_mach10", 303],
                   ["xi_0.5_mach100", 300],
                   ["xi_0.5_mach11.3", 300],
                   ["xi_0.5_mach160", 301],
                   ["xi_0.5_mach20", 301],
                   ["xi_0.5_mach24.51", 300],
                   ["xi_0.5_mach2.5", 300],
                   ["xi_0.5_mach320", 301],
                   ["xi_0.5_mach3.4", 300],
                   ["xi_0.5_mach40", 300],
                   ["xi_0.5_mach5", 313],
                   ["xi_0.5_mach50", 300],
                   ["xi_0.5_mach70", 496],
                   ["xi_0.5_mach7.71", 300],
                   ["xi_0.5_mach80", 301],
                   ["xi_0.75_mach10.7", 300],
                   ["xi_0.75_mach22", 300],
                   ["xi_0.75_mach3.5", 300],
                   ["xi_0.75_mach45", 300],
                   ["xi_0.75_mach7.2", 300],
                   ["xi_0.75_mach97", 300],
                   ["xi_0_mach10", 300],
                   ["xi_0_mach100", 117],
                   ["xi_0_mach124", 300],
                   ["xi_0_mach15.6", 300],
                   ["xi_0_mach35.3", 300],
                   ["xi_0_mach5.5", 300],
                   ["xi_0_mach65", 313],
                   ["xi_1_mach21.78", 300],
                   ["xi_1_mach3.4", 300],
                   ["xi_1_mach43", 300],
                   ["xi_1_mach7", 300],
                   ["xi_1_mach80", 38],
                   ["xi_1_mach94", 300]]
rrr = re.compile(r'xi_(.*)_mach(.*)')
sim_dir_base = "/anvil/scratch/x-ux454321/Paper83/256_mach_grid"
product_dir_base = "/anvil/scratch/x-ux454321/Paper83/256_mach_grid/products"
full_list=[]
for sim, nframe in list_and_frames:
    full_list.append(sim)
    match = rrr.match(sim)
    #short = match.group(1)
    ms = float(match.group(2))
    ma = 0
    simulation.sim(sim, data_location="%s/%s"%(sim_dir_base,sim),
                   product_location="%s/%s"%(product_dir_base,sim),
                   ms=ms,ma=ma,framelist=list(range(1,nframe+1)))

