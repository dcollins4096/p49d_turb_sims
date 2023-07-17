from  GL import *

import simulation
import simulation_info.all_sims

import data_scrub_2.compute_all_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
import data_scrub_2.image_frb as imf
import plots.P3_all_spectra as p3
import plots.P1_plot_quan as p1
#this_simname="1_1"
#this_simname="4s_dev"
#this_simname="4s_dave"
this_list=['4_half','4_2']

if 1:
    caq.comp_all(this_list)
if 0:
    m2d.make_all(this_list)
if 0:
    maf.make_all(this_list)
if 0:
    m3d.make_spec(sim_colors.simlist)

if 0:
    p1.plot_quan([this_simname])
if 0:
    for sim in this_list:
        p1.plot_quan([sim])
if 0:
    imf.image(this_simname)
if 0:
    p3.plot_all_spectra(sim_colors.simlist)


