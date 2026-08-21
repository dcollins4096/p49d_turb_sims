from GL import *


import simulation_info.suite_liltest
reload(simulation_info.suite_liltest)
import simulation_info.suite_1
reload(simulation_info.suite_1)
import simulation_info.suite_2
reload(simulation_info.suite_2)
import simulation_info.suite_2_256
reload(simulation_info.suite_2_256)
import simulation_info.suite_p83
import simulation_info.suite_3
reload(simulation_info.suite_3)
#import simulation_info.suite_4 as suite_4
#import simulation_info.suite_5 as suite_5
import simulation_info.suite_6_256 as suite_6
import simulation_info.suite_7_128 as suite_7
import simulation_info.suite_8_256 as suite_8
import simulation_info.brano_suite as suite_brano
import simulation_info.mach_grid as mach_grid
import simulation_info.brano_suite as brano_suite


lists={}
lists['brano'] = suite_brano.full_list
lists['suite1']=simulation_info.suite_1.simlist
lists['suite1a'] = ['half_half', 'half_1', 'half_2', 
                    '1_half', '1_1', '1_2', 
                    '2_half', '2_1', '2_2',
                    '3_half', '3_1', '3_2']
lists['suite1b'] = ['4_half', '4_1', '4_2', 
                    '5_half', '5_1', '5_2', 
                    '6_half', '6_1', '6_2']
lists['p83'] = ['4.7_3','4.7_4','x27_Ms8.0_Ma0.0_256']

lists['suite2'] = simulation_info.suite_2.long_simlist
lists['suite3'] = simulation_info.suite_3.long_simlist
lists['p83'] = ['4.7_3','4.7_4', '4_3','4_4', '2_3','2_4', '3_3','3_4','5_3','5_4']
lists['small'] = ['run2','run3']
#lists['suite4'] = [suite_4.long_from_key['d%02d'%n] for n in range(32)]
#lists['suite5'] = [suite_5.long_from_key['d%02d'%n] for n in range(16)]
lists['suite6'] = suite_6.longnamelist
lists['suite7'] = suite_7.longnamelist
lists['mach_grid'] = mach_grid.full_list
lists['brano_suite'] = brano_suite.full_list
