from GL import *
import simulation
reload(simulation)
import sim_colors

sim='4.7_3'
data = '/data/cb1/Projects/P83_MHDTurb/b03_Ms4.7_Ma3_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/b03_Ms4.7_Ma3_512'
ms = 4.7
ma = 3
color = 'r'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='4.7_4'
data = '/data/cb1/Projects/P83_MHDTurb/c03_Ms4.7_Ma4_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/c03_Ms4.7_Ma4_512'
ms = 4.7
ma = 4
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)
