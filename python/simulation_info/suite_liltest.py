from GL import *
import simulation
reload(simulation)


data_root = "/scratch/00369/tg456484/Paper49"
product_root = "/scratch/00369/tg456484/Paper49/Products"

simulation.sim("4s_dave", data_location=data_root+"/FourSmall_dave_vs", product_location=product_root+"/4s_dave", ms=4, ma=1,color='r',linestyle='--',marker='+',
               framelist=range(46))
               #framelist=range(0,46))
simulation.sim("4s_dev", data_location=data_root+"/FourSmall_enzo-dev", product_location=product_root+"/4s_dev", ms=4, ma=1,color='r',linestyle='--',marker='+',
               framelist=range(52))
               #framelist=range(0,46))
