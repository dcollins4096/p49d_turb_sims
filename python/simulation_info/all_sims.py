from GL import *


import simulation_info.suite_liltest
reload(simulation_info.suite_liltest)
import simulation_info.suite_1
reload(simulation_info.suite_1)

lists={}
lists['suite1']=simulation_info.suite_1.simlist
lists['suite1a'] = ['half_half', 'half_1', 'half_2', 
                    '1_half', '1_1', '1_2', 
                    '2_half', '2_1', '2_2',
                    '3_half', '3_1', '3_2']
lists['suite1b'] = ['4_half', '4_1', '4_2', 
                    '5_half', '5_1', '5_2', 
                    '6_half', '6_1', '6_2']

