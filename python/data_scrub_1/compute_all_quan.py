print('start')
from GL import *

import sim_colors

import compute_avg_quantities as comp_avg
reload(comp_avg)

#base_directory = sim_colors.cloudbreak_base
#output_directory_base = base_directory + "/Products/"
base_directory = sim_colors.stampede_run_base
output_directory_base = base_directory + "/Products/"
if len(sys.argv) > 1:
    simlist = [sys.argv[-1]]
else:
    simlist=sim_colors.simlist
print('simlist',simlist)
for sim in simlist:
    print("QUAN sim",sim)
    frames = sim_colors.framedict[sim]
    for frame in frames:
        comp_avg.make_quan(base_directory+sim,frame,out_directory=output_directory_base+sim )

