from dtools.starter1 import *

import shot
reload(shot)
import equal_probability_binner as ep
import physical_values as phys
reload(phys)

names = ['r60']#,'r120', 'r0']
if 'devices' not in dir() or True:
    devices={}
if 'r0' in  devices:
    del devices['r0']
bump_range={'r0':[175,300], 'r60':[175,300], 'r120':[175,300]}
vel_cut = {'r120':[325,450], 'r60':[325,450], 'r0':[410,550]}
#shock points should surround the foot of the shock and go past the peak.
shock_points = {'r60':[325,375, 650], 'r0':[325,375,650], 'r120':[325,375,650]}
shock_cut = {'r0':[550,750], 'r60':[475,675], 'r120':[475,675]}
for name in names:
    if name not in devices:
        tmp=shot.device(name, lines=[200,400], model=2,smooth=3)
        tmp.image_density(fname = 'image_shot_%s'%name)
        break
        x_off=None
        if name == 'r0':
            shift_60 = devices['r60'].shift_x
            shift_120 = devices['r120'].shift_x
            x_off = 0.5*(shift_60+shift_120)
        tmp.bumper(bump_range[name],fix_shift_x=x_off)#,fname='bumper_%s.pdf'%name)
        tmp.compute_velocity( vel_cut[name])#,fname = 'velocity_%s'%name)
        tmp.sigma_rho(shock_cut[name])#, fname = 'sigma_rho_%s'%name)
        tmp.csound(mean_density=means[name], fname='csound_%s'%name)
        #tmp.sigma_v( mean_density=means[name], fname = 'atwood_%s'%name)
        devices[name]=tmp
