from dtools.starter1 import *
import regions
reload(regions)
import physical_values as phys
reload(phys)

import horizontal_distance as horz
reload(horz)


import shot
reload(shot)

d_120=shot.device('r120')
rng=[175,300]
dx, drho = d_120.bumper(rng)

