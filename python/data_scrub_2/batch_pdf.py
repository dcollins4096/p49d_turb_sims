

from  GL import *

import simulation
import simulation_info.all_sims as all_sims

import data_scrub_2.make_all_pdf as mapdf
reload(mapdf)
sim_list=all_sims.lists['suite1']
mapdf.make_pdf(sim_list)
