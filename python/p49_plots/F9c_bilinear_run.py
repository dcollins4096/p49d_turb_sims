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
import itertools

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
    field_list = ['avg_clee','avg_clbb','avg_cltt','avg_v','avg_d','avg_h']
    #field_list = ['avg_d']
    #field_list = ['avg_clee','avg_clbb','avg_cltt']
    for nf,field in enumerate(field_list):
        msarr=[]
        maarr=[]
        amparr=[]
        slopearr=[]
        prob=False

        for sim in sim_colors.simlist:
            #if sim.startswith('5'):
            #    print('OOP',sim)
            #    continue
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
        herd[field + "s"] = bilinear.beefitter(field, msarr, maarr, slopearr)
        herd[field + "a"] = bilinear.beefitter(field, msarr, maarr, np.log(amparr))

field_order=['avg_clees', 'avg_cleea', 'avg_clbbs', 'avg_clbba', 'avg_cltts', 'avg_cltta', 'avg_vs', 'avg_va', 'avg_ds', 'avg_da', 'avg_hs', 'avg_ha']
field_order = ['avg_ds','avg_vs','avg_hs','avg_cltts','avg_clees','avg_clbbs',
               'avg_da','avg_va','avg_ha','avg_cltta','avg_cleea','avg_clbba']
bilinear.write_tex( herd, '%s/table1.tex'%plotdir, field_order)

#herd['avg_d'].plot2()
if 0:
    #predict
    truth = {'avg_clee':-2.4, 'avg_clbb':-2.5, 'avg_cltt':-2.6}
    truth_fields=truth.keys()
    if 0:
        Ntake = -1
        print('Test with %d'%Ntake)
        for A, B in itertools.combinations(truth_fields,2):
            t1 = herd[A].Q[Ntake]
            t2 = herd[B].Q[Ntake]
            Ms, Ma = bilinear.pairwise(herd[A], herd[B], t1,t2)
            print(Ms, Ma)
        print( msarr[Ntake], maarr[Ntake])

    print('real values')
    pred_ms = []; pred_ma = []
    for A, B in itertools.combinations(truth_fields,2):
        Ms, Ma = bilinear.pairwise(herd[A], herd[B], truth[A],truth[B])
        print(Ms,Ma)
        pred_ms.append(Ms)
        pred_ma.append(Ma)


    fig,ax=plt.subplots(1,1)
    for ms,ma in zip(pred_ms,pred_ma):
        ax.scatter(ms,ma, c='k')
    mean_ms,mean_ma=np.mean(pred_ms), np.mean(pred_ma)
    ax.scatter(mean_ms,mean_ma, c='k',marker='*')
    title=''
    for q in truth:
        title += '%s %0.2f '%(q[-2:],truth[q])
    title += 'Ms %0.2f Ma %0.2f'%(mean_ms,mean_ma)
    print(title)
    ax.set(xlabel=r'$M_s$', ylabel=r'$M_A$',title=title)
    fig.tight_layout()
    fig.savefig('%s/linear_prediction.pdf'%plotdir)

