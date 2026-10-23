from  GL import *

import simulation
import simulation_info
import simulation_info.all_sims as all_sims
import simulation_info.suite_5 as s5
import simulation_info.suite_7_128 as s7

import data_scrub_2.make_velocity_centroid as mvc

this_list = [simulation_info.suite_7_128.longname_from_key[sys.argv[1]]]
if 1:
    mvc.image_centroid(this_list)
if 0:
    m3d.make_spec(this_list)

#import plots.P2_image_all as p2
#for sim in this_list:
#    p2.image(sim)
