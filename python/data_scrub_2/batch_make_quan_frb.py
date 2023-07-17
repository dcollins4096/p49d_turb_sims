from  GL import *

#import simulation_info.suite_liltest
#reload(simulation_info.suite_liltest)
#import simulation_info.suite_1
#reload(simulation_info.suite_1)
import simulation
import simulation_info.all_sims as all_sims

import data_scrub_2.make_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf

if len(sys.argv) == 1:
    print("Please select 4, 5, 6")
    sys.exit(0)
#this_simname="1_1"
#this_simname="4s_dev"
#this_simname="4s_dave"
if sys.argv[1] == "4":
    this_list=['4_half','4_1','4_2']
elif sys.argv[1] == "5":
    this_list=['5_half','5_1','5_2']
elif sys.argv[1] == "6":
    this_list=['6_half','6_1','6_2']
elif sys.argv[1] == 'suite1':
    this_list = all_sims.lists['suite1']
else:
    print("Pick One.")
    sys.exit(0)
if 1:
    caq.comp_all(this_list)
if 1:
    maf.make_all(this_list)
if 0:
    m3d.make_spec(this_list)
