from  GL import *

#import simulation_info.suite_liltest
#reload(simulation_info.suite_liltest)
#import simulation_info.suite_1
#reload(simulation_info.suite_1)
import simulation
import simulation_info.suite_2
import simulation_info.all_sims as all_sims

import data_scrub_2.make_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
#this_simname="1_1"
#this_simname="4s_dev"
#this_simname="4s_dave"
#this_list = simulation_info.suite_2.list_from_key[sys.argv[1]]
this_list = [sys.argv[1]]

if 1:
    m3d.make_spec(this_list)
if 0:
    m2d.make_all(this_list)

import plots.P3_all_spectra as p3
if 1:
    #frames can be "all" or "ann"
    p3.plot_all_spectra(this_list, all_or_ann='ann')
