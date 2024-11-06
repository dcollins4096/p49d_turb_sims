from  GL import *

import simulation
import simulation_info.all_sims as all_sims
import simulation_info.suite_2

import data_scrub_2.make_cubes as mc

if len(sys.argv) == 1:
    print("Please select 4, 5, 6")
    sys.exit(0)

this_list = simulation_info.suite_2.list_from_key[sys.argv[1]]

mc.make_cubes(this_list)
