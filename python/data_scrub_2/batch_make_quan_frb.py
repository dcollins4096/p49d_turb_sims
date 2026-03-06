from  GL import *

import simulation
import simulation_info.all_sims as all_sims

import data_scrub_2.make_quan as caq
reload(caq)
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
#this_simname="1_1"
#this_simname="4s_dev"
#this_simname="4s_dave"
if 0:

    if len(sys.argv) == 1:
        print("Please select 4, 5, 6")
        sys.exit(0)
    if sys.argv[1] == "4":
        this_list=['4_half','4_1','4_2']
    elif sys.argv[1] == "5":
        this_list=['5_half','5_1','5_2']
    elif sys.argv[1] == "6":
        this_list=['6_half','6_1','6_2']
    elif sys.argv[1] == 'suite1':
        this_list = all_sims.lists['suite1']
    elif sys.argv[1] == '1':
        this_list = ['1_half','1_1','1_2']
    elif sys.argv[1] == '2':
        this_list = ['2_half','2_1','2_2']
    elif sys.argv[1] == '3':
        this_list = ['3_half','3_1','3_2']
    elif sys.argv[1] == 'half':
        this_list = ['half_half','half_1','half_2']
    else:
        print("Pick One.")
        sys.exit(0)

#this_list = ['half_half']
#this_list = [sys.argv[1]]
#this_list = all_sims.lists['suite1'][::-1]
#this_list = ['run3']
#this_list = ['6_2']
this_list = all_sims.lists['suite4'][3::4]
print(this_list)

if 0:
    caq.comp_all(this_list)
if 0:
    caq.energy_cleaner(this_list)
if 0:
    caq.comp_bulk(this_list)
if 0:
    caq.comp_Edot(this_list)
if 0:
    caq.comp_Ekin(this_list)

