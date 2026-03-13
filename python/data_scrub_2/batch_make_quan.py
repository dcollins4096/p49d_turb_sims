from  GL import *

import simulation
import simulation_info.all_sims as all_sims
import simulation_info.suite_4 as s4

import data_scrub_2.make_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
import simulation_info.suite_5
import simulation_info.suite_6_256
import simulation_info.suite_7_128
#this_list = simulation_info.suite_6_256.list_from_key[sys.argv[1]]
this_list = [sys.argv[1]]

#this_simname="1_1"
#this_simname="4s_dev"
#this_simname="4s_dave"
#if sys.argv[1] == 'ba':
#    this_list = ['ba_Ms2.0_Ma0.5_256']
#elif sys.argv[1] == 'aa':
#    this_list = ['aa_Ms2.0_Ma0.5_512']
#else:
#    print("Pick One.")
#    sys.exit(0)
if 1:
    caq.comp_all(this_list)
if 0:
    maf.make_all(this_list)
if 0:
    m3d.make_spec(this_list)

import plots.P1_plot_quan as p1
p1.plot_quan(this_list)
