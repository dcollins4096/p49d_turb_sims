
from GL import *
import queb3
reload(queb3)
import davetools as dt
import sim_colors
reload(sim_colors)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)
from scipy.optimize import curve_fit

#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
spectra_dict = rs.spectra_dict
import read_avg_quan as raq
quan3=raq.quan3

import bilinear
reload(bilinear)

if 1:
    Line1 = [1,2.1,3.4]
    Line2 = [1.5,1.5,1.5]
    Ms = [0.5,1,2]
    Ma = [0.5,1,2]
    def fun( line, ms, ma):
        return line[0]+line[1]*ms+line[2]*ma

    Q1 =nar([ fun(Line1, ms, ma) for ma in Ma for ms in Ms])
    Q2 =nar([ fun(Line2, ms, ma) for ma in Ma for ms in Ms])
    ms_all =nar([ms for ma in Ma for ms in Ms])
    ma_all =nar([ma for ma in Ma for ms in Ms])

    bf1 = bilinear.beefitter( 'L1', ms_all, ma_all, Q1)
    bf2 = bilinear.beefitter( 'L2', ms_all, ma_all, Q2)

    t1 = fun(Line1, 0.62, 2.3)
    t2 = fun(Line2, 0.62, 2.3)
    found_s, found_a = bilinear.pairwise( bf1, bf2, t1, t2)
    print(found_s,found_a)
