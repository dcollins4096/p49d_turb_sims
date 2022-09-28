
from GL import *

import sim_colors

import compute_avg_quantities as comp_avg
reload(comp_avg)

base_directory = sim_colors.cloudbreak_base
output_directory_base = base_directory + "/Products/"
for sim in sim_colors.simlist:
    print("QUAN sim",sim)
    frames = sim_colors.framedict[sim]
    for frame in frames:
        comp_avg.make_quan(base_directory+sim,frame,out_directory=output_directory_base+sim )

