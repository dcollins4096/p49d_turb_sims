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

plotdir =  "/home/dccollins/PigPen"

#save bi-linear fit parameters to an h5 file
#function recording values for yaxis (perp to b-field)
if 1:
    herd = {}
    #field_list = ['avg_clee','avg_clbb','avg_cltt','avg_v','avg_d','avg_h']
    field_list = ['avg_clee','avg_clbb','avg_cltt']
    truth = {'avg_clee':-2.4, 'avg_clbb':-2.5, 'avg_cltt':-2.3}
    for nf,field in enumerate(field_list):
        msarr=[]
        maarr=[]
        amparr=[]
        slopearr=[]
        prob=False

        for sim in sim_colors.simlist:
            if (spectra_dict['y'][sim].slopes[field]==0 or spectra_dict['z'][sim].slopes[field] == 0):
                continue
            #msarr=np.append(msarr,sim_colors.Mst[sim])
            #maarr=np.append(maarr,sim_colors.Mat[sim])
            #amparr=np.append(amparr,spectra_dict['y'][sim].amps[field])
            #slopearr=np.append(slopearr,spectra_dict['y'][sim].slopes[field])
            msarr=np.append(msarr,raq.quan3[sim]['msavg'])
            maarr=np.append(maarr,raq.quan3[sim]['maavg'])
            amparr=np.append(amparr,spectra_dict['z'][sim].amps[field])
            slopearr=np.append(slopearr,spectra_dict['z'][sim].slopes[field])
            #prob=True
        if (prob):
            do_morebifit("%s_slopes"%field,msarr,maarr,slopearr)
            do_morebifit("%s_amps"%field,msarr,maarr,amparr)
        herd[field] = bilinear.beefitter(field, msarr, maarr, slopearr)


    Ntake = -1
    print('Test with %d'%Ntake)
    for A, B in itertools.combinations(field_list,2):
        t1 = herd[A].Q[Ntake]
        t2 = herd[B].Q[Ntake]
        Ms, Ma = bilinear.pairwise(herd[A], herd[B], t1,t2)
        print(Ms, Ma)
    print( msarr[Ntake], maarr[Ntake])

    print('real values')
    truth = {'avg_clee':-2.4, 'avg_clbb':-2.5, 'avg_cltt':-2.7}
    for A, B in itertools.combinations(field_list,2):
        Ms, Ma = bilinear.pairwise(herd[A], herd[B], truth[A],truth[B])
        print(Ms, Ma)


