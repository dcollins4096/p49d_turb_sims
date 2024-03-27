#
# Main python packages of my design.
#
from starter1 import *

path_list=['tools','OtherPython']
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]

import yt
import davetools as dt
import queb3
import sim_colors
reload(sim_colors)
reload(queb3)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)
import data_locations as dl
root4pi = np.sqrt(4*np.pi)
plotdir = "%s/PigPen"%os.environ['HOME']
