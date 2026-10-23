from  GL import *

import simulation
import simulation_info
import simulation_info.all_sims as all_sims

import data_scrub_2.make_ppv as ppv

this_list = [sys.argv[1]]

ppv.make_ppv(this_list)

