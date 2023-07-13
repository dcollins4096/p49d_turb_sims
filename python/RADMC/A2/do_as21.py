from GL import *
import A2 
reload(A2)
from A2 import problem_setup
repo_dir = "%s/repos/p49d_turb_sims/python/A2/"%os.environ['HOME']
output_dir="/data/cb1/Projects/P49_EE_BB/RMC_128/1_1"
output_dir="/data/cb1/Projects/P49_EE_BB/RMC_128/as21_M0.6_MA0.3_64"
maker = problem_setup.rmc_filemaker(repo_dir=repo_dir, output_dir=output_dir, dust_model='simple_wrong')
#maker.write_dust_parameters()
test_set_1 = "/data/cb1/Projects/P49_EE_BB/Simulations_128_downsample/1_1/DD0011/data0011"
test_set_2 = "/data/cb1/Projects/P49_EE_BB/as21_M0.6_MA0.3_64/DD1080/data1080"
maker.write_radmc_from_enzo(test_set_2)


