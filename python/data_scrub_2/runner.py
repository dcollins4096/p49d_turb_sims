from  GL import *

import simulation_info.suite_liltest
reload(simulation_info.suite_liltest)
import simulation_info.suite_1
reload(simulation_info.suite_1)

import data_scrub_2.compute_all_quan as caq
#reload(caq)
#this_simname="1_1"
#this_simname="4s_dev"
this_simname="4s_dave"
if 0:
    caq.comp_all(this_simname)
import plots.P1_plot_quan as p1
reload(p1)
if 0:
    p1.plot_quan([this_simname])
if 0:
    for sim in sim_colors.simlist:
        p1.plot_quan([sim])

import data_scrub_2.make_all_frb as maf
reload(maf)
if 0:
    maf.make_all([this_simname])
import data_scrub_2.image_frb as imf
reload(imf)
if 1:
    imf.image(this_simname)


