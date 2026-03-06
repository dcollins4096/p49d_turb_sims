
from GL import *
import simulation
reload(simulation)
import sim_colors
import re

full_list=['d00_Ms1.0_Ma0.5_512', 'd01_Ms1.0_Ma1.5_512', 'd02_Ms1.0_Ma3.0_512', 'd03_Ms1.0_Ma0.0_512', 'd04_Ms2.0_Ma0.5_512', 'd05_Ms2.0_Ma1.5_512', 'd06_Ms2.0_Ma3.0_512', 'd07_Ms2.0_Ma0.0_512', 'd08_Ms3.0_Ma0.5_512', 'd09_Ms3.0_Ma1.5_512', 'd10_Ms3.0_Ma3.0_512', 'd11_Ms3.0_Ma0.0_512', 'd12_Ms4.0_Ma0.5_512', 'd13_Ms4.0_Ma1.5_512', 'd14_Ms4.0_Ma3.0_512', 'd15_Ms4.0_Ma0.0_512', 'd16_Ms5.0_Ma0.5_512', 'd17_Ms5.0_Ma1.5_512', 'd18_Ms5.0_Ma3.0_512', 'd19_Ms5.0_Ma0.0_512', 'd20_Ms6.0_Ma0.5_512', 'd21_Ms6.0_Ma1.5_512', 'd22_Ms6.0_Ma3.0_512', 'd23_Ms6.0_Ma0.0_512', 'd24_Ms7.0_Ma0.5_512', 'd25_Ms7.0_Ma1.5_512', 'd26_Ms7.0_Ma3.0_512', 'd27_Ms7.0_Ma0.0_512', 'd28_Ms8.0_Ma0.5_512', 'd29_Ms8.0_Ma1.5_512', 'd30_Ms8.0_Ma3.0_512', 'd31_Ms8.0_Ma0.0_512']

rrr = re.compile(r'(...)_Ms(...)_Ma(...)_512')
for sim in full_list:
    match = rrr.match(sim)
    print(match.group(3))
    name = match.group(0)
    short = match.group(1)
    ms = match.group(2)
    ma = match.group(3)
    simulation.sim(short, data_location="%s/%s"%(dl.sim_dir_base,name),
                   product_location="%s/%s"%(dl.product_dir_base,name),
                   ms=ms,ma=ma,framelist=list(range(102)))

