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
from collections import defaultdict
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
import read_stuff as rs
def cube(x,a,b,c,d):
    return a + b*x + c*x**2 +d*x**3
def cube_prime(x,a,b,c,d):
    return  b + 2*c*x +3*d*x**2
def cube_prime_prime(x,a,b,c,d):
    return   2*c +6*d*x
def quad(x,a,b,c):
    return a + b*x + c*x**2
def line(x,a,b):
    return a + b*x
def minmax(alpha,x,y):
    v = x**alpha*y
    return nar([min(v),max(v)])
def broken_powerlaw(x, x0, A, m1, m2):
    output=np.zeros_like(x)
    ok = x-x0>0
    output[  ok ] = (x-x0)[ok]*m1 + A
    output[ ~ok ] = (x-x0)[~ok]*m2 + A
    return output
                        
class slope_bucket():
    def __init__(self):
        self.slope_array=[]

        self.cube_slope=[]
        self.cube_amp=[]
        self.cube_slope_list=[]
        self.line_slope_list=[]
        self.cube_amp_list=[]
        self.cube_param_list=None
        self.k_flat_list=[]
        self.ok_list=[]
        self.chi_l_list=[]
        self.chi_c_list=[]

        self.line_list=[]
        self.cube_list=[]
        self.xvals=[]

        self.AlphaF=None
        self.AlphaPeak=None

        self.broken_spec=None
        self.broken_xval=None
        self.broken_fits=None



field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
spectra_dict=rs.spectra_dict
if 1:
    #fill objects with slope array
    for nf,field in enumerate(field_list):
        for ns,sim in enumerate(sim_colors.simlist):
            print(field,sim)
            for LOS in 'xyz':
                spectra_dict[LOS][sim].slopes3[field]=slope_bucket()
                mah_bucket=spectra_dict[LOS][sim].slopes3[field]
                proj=rs.proj_dict[LOS][sim]
                x0_list = proj.lcent[5:10:2]
                x1_list = proj.lcent[14:51:2]



                max_npoints=np.argmax(x1_list)-np.argmin(x0_list)+1
                mah_bucket.x0_list=x0_list
                mah_bucket.x1_list=x1_list
                specf_raw = spectra_dict[LOS][sim].spectra[field][5:]
                specf = 10**gaussian_filter1d( np.log10(specf_raw), 1, mode='nearest')
                mah_bucket.specf = specf
                mah_bucket.specf_raw = specf_raw
                mah_bucket.lcentf = proj.lcent[5:]
                lcentf=mah_bucket.lcentf

                mah_bucket.broken_spec = specf
                mah_bucket.broken_xval = lcentf
                fit, cov = curve_fit( broken_powerlaw, np.log10(lcentf), np.log10(specf))
                mah_bucket.broken_fits=fit

                for n0,x0 in enumerate(mah_bucket.x0_list):
                    for n1,x1 in enumerate(mah_bucket.x1_list):

                        ok = (lcentf>=x0)*(lcentf<=x1)
                        xvals = lcentf[ok]
                        spec=specf[ok]
                            
                        popt_1, pcov_1= curve_fit( line, np.log10(xvals), np.log10(spec))
                        popt_3, pcov_3= curve_fit( cube, np.log10(xvals), np.log10(spec))
                        a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                        slope = b3 - c3**2/(3*d3)
                        value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                        mah_bucket.line_slope_list.append(popt_1)
                        mah_bucket.cube_slope_list.append(slope)
                        mah_bucket.cube_amp_list.append(value)
                        if  mah_bucket.cube_param_list is None:
                            mah_bucket.cube_param_list=popt_3
                        else:
                            mah_bucket.cube_param_list=np.vstack([mah_bucket.cube_param_list,popt_3])
                        mah_bucket.ok_list.append(ok)

                        this_c=10**cube(np.log10(xvals), popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                        this_l=10**line(np.log10(xvals), popt_1[0],popt_1[1])
                        mah_bucket.xvals.append(xvals)
                        mah_bucket.line_list.append(this_l)
                        mah_bucket.cube_list.append(this_c)
#
                        chi_c = ((np.log10(this_c)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                        chi_l = ((np.log10(this_l)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                        mah_bucket.chi_l_list.append(chi_l)
                        mah_bucket.chi_c_list.append(chi_c)

                #bins = np.mgrid[-5:0:16j]
                mah_bucket.cube_param_list=nar(mah_bucket.cube_param_list)
                bins = np.mgrid[-5:0:51j]
                hist, bins = np.histogram(mah_bucket.cube_slope_list, bins=bins)
                most_probable_bin = np.argmax(hist)
                bL = bins[most_probable_bin]
                bR = bins[most_probable_bin+1]
                SL = nar(mah_bucket.cube_slope_list)
                AL = nar(mah_bucket.cube_amp_list)
                AlphaC = SL[ (SL >= bL)*(SL <= bR)].mean()
                AmpC = AL[ (SL >= bL)*(SL <= bR)].mean()
                mah_bucket.AlphaC = AlphaC
                mah_bucket.AmpC = AmpC

                #do it again for the linear slope
                bins = np.mgrid[-5:0:51j]
                SL = nar(nar(mah_bucket.line_slope_list)[:,1])
                hist1, bins1 = np.histogram(SL, bins=bins)
                most_probable_bin = np.argmax(hist1)
                bL = bins1[most_probable_bin]
                bR = bins1[most_probable_bin+1]
                most_probable =  (SL >= bL)*(SL <= bR)

                mah_bucket.AlphaPeak = SL[most_probable].mean()
