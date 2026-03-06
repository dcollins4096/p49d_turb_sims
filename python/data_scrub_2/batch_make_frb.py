from  GL import *

import simulation
import simulation_info
import simulation_info.all_sims as all_sims
#import simulation_info.suite_5 as s5
#import simulation_info.suite_7_128 as s7

import data_scrub_2.make_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
#this_list = simulation_info.suite_2.list_from_key[sys.argv[1]]
#this_list = [s5.long_from_key[sys.argv[1]]]

this_list = ['half_half']
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
    maf.make_all(this_list)
if 0:
    m3d.make_spec(this_list)

#import plots.P2_image_all as p2
#for sim in this_list:
#    p2.image(sim)
