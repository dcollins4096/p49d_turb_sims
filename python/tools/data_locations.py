from GL import *

cloudbreak_base = "/data/cb1/Projects/P49_EE_BB/"
stampede_run_base = "/scratch/00369/tg456484/Paper49/RUNNING/"
stampede_base = "/scratch/00369/tg456484/Paper49/RUNNING/"
cloudbreak_128 = "/data/cb1/Projects/P49_EE_BB/Downsample128"
p58_dir = "/data/cb1/Friends/P58_synchrotron"
plotdir="./plots_to_sort"
M
if os.environ['machine']=='stampede2':
    sim_dir_base = stampede_base
    product_dir_base = sim_dir_base + "/Products/"
else:
    sim_dir_base = cloudbreak_base
    product_dir_base = base_directory + "/Products/"
