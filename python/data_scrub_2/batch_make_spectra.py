from  GL import *

#import simulation_info.suite_liltest
#reload(simulation_info.suite_liltest)
#import simulation_info.suite_1
#reload(simulation_info.suite_1)
import simulation
import simulation_info.all_sims

import data_scrub_2.compute_all_quan as caq
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
import data_scrub_2.image_frb as imf
import plots.P3_all_spectra as p3
import plots.P1_plot_quan as p1
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
else:
    print("Pick 4 5 6")
    sys.exit(0)

if 1:
    m3d.make_spec(this_list)

