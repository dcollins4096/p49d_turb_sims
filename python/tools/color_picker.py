
from GL import *
import sim_colors
import davetools as dt

rm = dt.rainbow_map(30)
plt.clf()
for n in range(30):
    plt.scatter(n,n,c=[rm(n)])
plt.savefig('/home/dccollins/PigPen/colortest.png')

