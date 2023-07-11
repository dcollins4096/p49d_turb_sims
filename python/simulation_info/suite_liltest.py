from GL import *
import simulation
reload(simulation)


data_root = "/scratch/00369/tg456484/Paper49/FourSmall_dave_vs"
product_root = "/scratch/00369/tg456484/Paper49/Products"

simulation.sim("4s_dave", data_location=data_root, product_location=product_root+"/4s_dave", ms=4, ma=1,color='r',linestyle='--',marker='+',
               framelist=range(46))
               #framelist=range(0,46))
