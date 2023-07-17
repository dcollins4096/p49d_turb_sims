from GL import *

import simulation
import simulation_info.all_sims as all_sims

import plots.P1_plot_quan as p1
import plots.P2_image_all as p2
import plots.P3_all_spectra as p3
import plots.P4_spectra_time as p4
reload(p1)
reload(p2)
reload(p3)
reload(p4)

sim_list = all_sims.lists['suite1b']


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
    #frames can be "all" or "ann"
    p3.plot_all_spectra(sim_list, all_or_ann='ann')

if 0:
    p4.slope_time(sim_list)
