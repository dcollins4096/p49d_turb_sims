from starter1 import *
import queb3
import sim_colors
reload(sim_colors)
reload(queb3)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)

path_list=['tools','OtherPython']
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]

