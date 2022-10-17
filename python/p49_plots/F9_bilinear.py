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

plotdir =  "/home/dccollins/PigPen"

def bifit(x, a, b, c, d):
    #x0 is Ms x1 is Ma
    return a+ b*x[0] + c*x[1] + d*x[0]*x[1]
def do_bifit(xdata, ydata, p0):
    fitParams, fitCovariances = curve_fit(bifit, xdata, ydata, p0)
    print(' fit coefficients:\n', fitParams)
    return fitParams, fitCovariances
def do_morebifit(fieldname, ms, ma, myval):
    x=np.row_stack([ms.flatten(),ma.flatten()])
    y=myval.flatten()
    p0 = [0,1,1,0]
    print(fieldname)
    fitParams, fitCovariances= do_bifit(x,y, p0)
    f = fitParams
    cov = fitCovariances
#    ys = f[0] + f[1]*myx + f[2]*myy+f[3]*myx*myy
    Fptr = h5py.File("bi_fit_params.h5",'a')
    Fptr.create_dataset(fieldname, data=f)
    Fptr.close()
    Fptr = h5py.File("bi_fit_covariances.h5",'a')
    Fptr.create_dataset(fieldname+"_cov", data=cov)
    Fptr.close()

#save bi-linear fit parameters to an h5 file
#function recording values for yaxis (perp to b-field)
if 1:
    for nf,field in enumerate(['avg_clee','avg_clbb','avg_cltt','avg_v','avg_d','avg_h']):
        msarr=[]
        maarr=[]
        amparr=[]
        slopearr=[]
        prob=False
        fptr = h5py.File('bi_fit_params.h5','w')
        fptr.close()
        fptr = h5py.File('bi_fit_covariances.h5','w')
        fptr.close()

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
            prob=True
        if (prob):
            do_morebifit("%s_slopes"%field,msarr,maarr,slopearr)
            do_morebifit("%s_amps"%field,msarr,maarr,amparr)

