
from GL import *

import sim_colors

import compute_avg_quantities as comp_avg
reload(comp_avg)

for sim in ['5_half','5_1','5_2','5_3']:
    print("QUAN",sim)
    frames = sim_colors.framedict[sim]
    for frame in frames:
        comp_avg.make_quan(sim_colors.cloudbreak_base+sim,frame)
    break

