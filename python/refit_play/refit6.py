"""
refit 5.  
Fits are done in refit3.
Plotted here are 
* left: compensated spectra; smoothed spectra; Fixed Fit line
* right: Distribution of alpha
"""


from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(sim_colors)
reload(queb3)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)

field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
#simlist = ['half_half','1_half','2_half','3_half']#['1_half']#['half_2', 'half_half']#['half_half']#['1_1']#['2_2'] #['half_half'] # ['half_half']#,'2_2']
#simlist = ['3_half', 'half_half', '3_2']
simlist=['3_1']
field_list = ['avg_clbb']
do_all_subplots=False


#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
#reload(read_stuff)  
#spectra_dict = rs.spectra_dict
quan3 = rs.quan3
import refit3
#reload(refit3)

spectra_dict = refit3.spectra_dict

from scipy.optimize import curve_fit

#
# Spectra plots
#

plotdir=os.environ['HOME']+"/PigPen"

from scipy.ndimage import gaussian_filter1d
def compensate_plot(ax,alpha,x,y,**kwargs):
    ax.plot( x, x**alpha*y, **kwargs)
def compensate_scatter(ax,alpha,x,y,**kwargs):
    ax.scatter( x, x**alpha*y, **kwargs)
if 1:
    #
    # Development code.  
    #
    #fill objects with slope array

    fig_all,ax_all=plt.subplots(1,1, figsize=(8,12))
    for nf,field in enumerate(field_list):
        LOS_LIST = 'y'
        if field in ['avg_d','avg_v','avg_h']:
            LOS_LIST='x'

        for ns,sim in enumerate(simlist):
            def broken_powerlaw(x, x0_drp, A, m1, m2):
                #return x*m1+A
                #x0=np.log10(0.5)
                output=np.zeros_like(x)
                ok = x-x0>0
                output[  ok ] = (x-x0)[ok]*m1 + A
                output[ ~ok ] = (x-x0)[~ok]*m2 + A
                return output
            
            for LOS in LOS_LIST:
                mah_bucket = spectra_dict[LOS][sim].slopes3[field]
                print("DO", field, sim, LOS)
                proj=rs.proj_dict[LOS][sim]
                fit_range=proj.determine_fit_range()
                x0=fit_range[0]
                x1=fit_range[1]
                lcentf = mah_bucket.lcentf
                specf =  mah_bucket.specf
                specf_raw = mah_bucket.specf_raw

                mah_bucket.broken_spec = specf[:100]
                mah_bucket.broken_xval = lcentf[:100]
                fit, cov = curve_fit( broken_powerlaw, np.log10(mah_bucket.broken_xval), np.log10(mah_bucket.broken_spec))
                mah_bucket.broken_fits=fit



                # fixed fit
                ok = (lcentf>=x0)*(lcentf<=x1)
                xvals = lcentf[ok]
                spec=specf[ok]
                popt_fixed, pcov_fixed= curve_fit( refit3.line, np.log10(xvals), np.log10(spec))
                AlphaF = popt_fixed[1]
                AlphaComp1 = AlphaF

                #AlphaComp1 = 0 #AlphaF

                verts=dt.extents( refit3.minmax( -AlphaComp1, lcentf[ok], spec))
                verts=dt.extents( refit3.minmax( -AlphaComp1, lcentf, specf_raw))

                fig_sub,ax=plt.subplots(1,1,figsize=(12,8))
                ax0=ax
                #ax_sub=ax_sub_arr[0]
                #ax_sub_1=ax_sub_arr[1]
                #ax_sub_2=ax_sub_arr[2]

                #
                # Spectra plots
                #


                this_l = 10**broken_powerlaw( np.log10(mah_bucket.broken_xval), *mah_bucket.broken_fits)
                verts( refit3.minmax( -AlphaComp1, mah_bucket.broken_xval, this_l))
                compensate_plot(ax0, -AlphaComp1, lcentf, specf, c=[0]*3)
                compensate_plot(ax0, -AlphaComp1, lcentf, specf_raw, c=[0.8]*3)
                compensate_plot(ax0, -AlphaComp1,  mah_bucket.broken_xval, this_l,c='g')


                #
                # save
                #

                xlim = [proj.lcent.min(),proj.lcent.max()]
                dt.axbonk(ax0,xscale='log',yscale='log', ylim=verts.minmax, xlim=xlim, xlabel=r'$k$',ylabel=r'$P*k^{%0.1f}$'%AlphaComp1)
                fig_sub.savefig('%s/%s_%s_%s_broken_powerlaw.png'%(plotdir,sim,field,LOS))
