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

sim='x27_Ms8.0_Ma0.0_256'
data = '/anvil/scratch/x-ux454321/p83_turbulence/Athena/maker/x27_Ms8.0_Ma0.0_256'
product = '/anvil/scratch/x-ux454321/p83_turbulence/Athena/Products/%s'%sim
ms = 8
ma = 0
color = 'g'; line='--';marker='*'
framelist = list(range(1,99))
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist, all_frames=framelist, code='Athena')

sim='4_3'
data = '/data/cb1/Projects/P83_MHDTurb/d03_Ms4_Ma3_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/d03_Ms4_Ma3_512'
ms = 4
ma = 3
color = 'r'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='4_4'
data = '/data/cb1/Projects/P83_MHDTurb/e03_Ms4_Ma4_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/e03_Ms4_Ma4_512'
ms = 4
ma = 4
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='2_3'
data = '/data/cb1/Projects/P83_MHDTurb/f03_Ms2_Ma3_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/f03_Ms2_Ma3_512'
ms = 2
ma = 3
color = 'r'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='2_4'
data = '/data/cb1/Projects/P83_MHDTurb/g03_Ms2_Ma4_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/g03_Ms2_Ma4_512'
ms = 2
ma = 4
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='3_3'
data = '/data/cb1/Projects/P83_MHDTurb/h03_Ms3.0_Ma3.0_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/h03_Ms3.0_Ma3.0_512'
ms = 3
ma = 3
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='3_4'
data = '/data/cb1/Projects/P83_MHDTurb/i03_Ms3.0_Ma4.0_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/i03_Ms3.0_Ma4.0_512'
ms = 3
ma = 4
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='5_3'
data = '/data/cb1/Projects/P83_MHDTurb/j03_Ms5.0_Ma3.0_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/j03_Ms5.0_Ma3.0_512'
ms = 5
ma = 3
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)

sim='5_4'
data = '/data/cb1/Projects/P83_MHDTurb/k03_Ms5.0_Ma4.0_512'
product = '/data/cb1/Projects/P83_MHDTurb/Products/k03_Ms5.0_Ma4.0_512'
ms = 5
ma = 4
color = 'g'; line='--';marker='*'
framelist = [100]
simulation.sim(sim,data_location=data,product_location=product,ms=ms,ma=ma,color=color,linestyle=line,marker=marker,framelist=framelist)
