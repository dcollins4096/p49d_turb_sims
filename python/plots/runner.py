from GL import *

import simulation
import simulation_info.all_sims as all_sims

import plots.P1_plot_quan as p1
import plots.P2_image_all as p2
import plots.P3_all_spectra as p3
reload(p1)
reload(p2)
reload(p3)

sim_list = all_sims.lists['suite1']
for n,s in enumerate(sim_list):
    print(n,s)
sim_list = sim_list[18:]


if 0:
    for sim in sim_list:
        this_sim = simulation.corral[sim]
        print(this_sim.ann_frames)

if 0:
    #all plots in one pannel
    p1.plot_quan(sim_list)
if 0:
    #quan plot, each sim
    for sim in sim_list:
        p1.plot_quan([sim])
if 0:
    #12 panel image for each frame
    for sim in sim_list:
        p2.image(sim)

if 1:
    p3.plot_all_spectra(sim_list)

