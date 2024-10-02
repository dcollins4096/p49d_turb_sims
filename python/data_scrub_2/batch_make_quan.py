from  GL import *

import simulation
import simulation_info.all_sims as all_sims
import simulation_info.suite_2

import data_scrub_2.make_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf

this_list = simulation_info.suite_2.list_from_key[sys.argv[1]]
if 1:
    caq.comp_all(this_list)
if 0:
    maf.make_all(this_list)
if 0:
    m3d.make_spec(this_list)

import plots.P1_plot_quan as p1
p1.plot_quan(this_list)
