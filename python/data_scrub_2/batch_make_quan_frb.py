from  GL import *

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
if 1:
    caq.comp_all(this_list)
if 0:
    maf.make_all(this_list)
if 0:
    m3d.make_spec(this_list)

