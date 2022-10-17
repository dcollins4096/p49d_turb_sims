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

def cube(x,a,b,c,d):
    return a + b*x + c*x**2 +d*x**3
class slopes():
    def __init__(self,name):
        self.name=name
        self.spectra={}
        self.std={}
        self.std_s={}
        self.slopes={}
        self.amps={}
        self.res={}
        self.pack = None
        self.slopes3={}
    def fit(self,lcent,fitrange):
        for field in ['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']:

            slope,amp,res=plfit(lcent,self.spectra[field],fitrange)
            self.amps[field]=amp
            self.slopes[field]=slope
            self.res[field]=res
def pline(self,ax,field,norm=False,label=False,**kwargs):
    """Plots the line we found onto *ax*.
    Also update the label of the line to reflect the value
    of the slope"""
    m=self.slopes[field]
    a=self.amps[field]
    ellfit=1    
    if norm:
        ellfit = np.sqrt(self.fit_range[0]*self.fit_range[1])
    y0 = a*(self.fit_range[0]/ellfit)**m
    y1 = a*(self.fit_range[1]/ellfit)**m
    if label:
        label=kwargs.get('label','')
        label += r' $%0.1f$'%m
        kwargs['label']=label
    ax.plot( self.fit_range,[y0,y1],**kwargs)

#
# Read spectra
#

spectra_dict={}
proj_dict={}
shortprefix='time'
for axes in ['x','y','z']:
    spectra_dict[axes]={}
    proj_dict[axes]={}
    for i, sim in enumerate(sim_colors.simlist):
        spectra_dict[axes][sim]=slopes(sim) 

        sim_dir = "/data/cb1/Projects/P49_EE_BB/%s"%sim
        product_dir = "/data/cb1/Projects/P49_EE_BB/Products/%s"%sim
        longprefix='%s_%s_%s'%(sim,shortprefix,axes)
        spectra_fname = "avg_spectra_%s_%s.h5"%(sim,axes)
        pack = queb3.simulation_package( directory=sim_dir,frames=sim_colors.frames[sim],prefix=longprefix,
                                        product_directory=product_dir)
        proj=pack.read_queb(frame=sim_colors.frames[sim][0],ax=axes,bin_style='dx1')
        fitrange=proj.determine_fit_range()  #something better
        proj_dict[axes][sim]=proj
        proj.fitrange=fitrange
        if os.path.exists(spectra_fname):
            try:
                fptr=h5py.File(spectra_fname,'r')
                for field in fptr:
                    spectra_dict[axes][sim].spectra[field]=fptr[field][()]

                    if field in ['avg_cltt','avg_clbb','avg_clee']:
                        std_field = "std_%s"%(field.split('_')[1])
                        spectra_dict[axes][sim].std[field]=fptr[std_field][()]
                        std_s_field = "std_s_%s"%(field.split('_')[1])
                        spectra_dict[axes][sim].std_s[field]=fptr[std_s_field][()]
                        #print('wut', axes, sim, field, spectra_dict[axes][sim].std_s[field])


            except:
                raise
            finally:
                fptr.close()
        spectra_dict[axes][sim].fit(proj.lcent,fitrange)
        spectra_dict[axes][sim].lcent = proj.lcent
        spectra_dict[axes][sim].fit_range = fitrange    

#print("WTF",spectra_dict['y'][sim].std_s['avg_cltt'])


#
# Average slopes
#

