from GL import *

if os.environ['machine']=='stampede2':
    stampede_run_base = "/scratch/00369/tg456484/Paper49/RUNNING/"
    stampede_base = "/scratch/00369/tg456484/Paper49/RUNNING/"
    sim_dir_base = stampede_base
    product_dir_base = sim_dir_base + "/Products/"
elif os.environ['machine']=='cloudbreak':
    cloudbreak_base = "/data/cb1/Projects/P49_EE_BB/"
    cloudbreak_128 = "/data/cb1/Projects/P49_EE_BB/Downsample128"
    p58_dir = "/data/cb1/Friends/P58_synchrotron"
    sim_dir_base = cloudbreak_base
    product_dir_base = sim_dir_base + "/Products/"
elif os.environ['machine']=='teahupoo':
    sim_dir_base = "/data/cb1/Projects/P49_EE_BB/"
    product_dir_base = "%s/Products/"%sim_dir_base
    plotdir = os.environ['HOME']+"/PigPen/"


    
