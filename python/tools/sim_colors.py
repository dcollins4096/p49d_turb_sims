#color is Alfven
#line is Sonic
#marker is Sonic
from GL import *

cloudbreak_base = "/data/cb1/Projects/P49_EE_BB/"
cloudbreak_128 = "/data/cb1/Projects/P49_EE_BB/Downsample128"

#these need to stay in this order.
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2","5_half","5_1","5_2"])
#simlist=nar(["5_half"]); print('kludge in simlist')
#simlist=nar(["1_1"])
color_list=nar(['r','g','b','r','g','b','r','g','b','r','g','b','r','g','b'])
line_list=nar(['-','-','-','-.','-.','-.','--','--','--',':',':',':', '-','-','-'])
glyph_list = nar([a+b for a,b in zip(color_list,line_list)])
marker_list = nar(['.','.','.','s','s','s','^','^','^','*','*','*','o','o','o'])


framedict={
    "half_half":range(11,31),"half_1":range(11,31),"half_2":range(11,31),"1_half":range(11,31),"1_1":range(11,31),"1_2":range(11,31),"2_half":range(65,85),"2_1":range(11,31),"2_2":range(11,31),"3_half":range(72,91),"3_1":range(56,75),"3_2":range(20,40),"5_half":range(3,36),"5_1":range(4,27),"5_2":range(4,45),"5_3":range(5,58)}
#framedict['5_half']=range(11,13); print('kludge in framedict')

sim_ms = nar(['half','1','2','3','5'])
sim_ma = nar(['half','1','2'])
framelist=[framedict[sim] for sim in simlist]

plot_order=[]
for ma in sim_ma:
    for ms in sim_ms:
        sim="%s_%s"%(ms,ma)
        plot_order.append(sim)

color=dict(zip(simlist,color_list))
linestyle = dict(zip(simlist,line_list))
marker = dict(zip(simlist,marker_list))
frames = dict(zip(simlist,framelist))
glyph = dict(zip(simlist,glyph_list))

def vals_from_sim(sim):
    ms,ma = sim.split("_")
    if ms == 'half':
        ms = 0.5
    if ma == 'half':
        ma = 0.5
    ms=float(ms)
    ma=float(ma)
    return ms,ma

ms_list=[]
ma_list=[]
for sim in simlist:
    ms,ma = vals_from_sim(sim)
    ms_list.append( ms)
    ma_list.append(ma)
ms_list=nar(ms_list)
ma_list=nar(ma_list)
Ms = dict(zip(simlist,ms_list))
Ma = dict(zip(simlist,ma_list))


