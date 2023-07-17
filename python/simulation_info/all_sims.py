from GL import *


import simulation_info.suite_liltest
reload(simulation_info.suite_liltest)
import simulation_info.suite_1
reload(simulation_info.suite_1)

lists={}
lists['suite1']=simulation_info.suite_1.simlist
