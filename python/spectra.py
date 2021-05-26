
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

#
# Read primitives
# 

class slopes():
    def __init__(self,name):
        self.name=name
        self.spectra={}
        self.slopes={}
        self.amps={}
        self.res={}
        self.pack = None
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
if 'quand' not in dir():
    quand = {}
for i, sim in enumerate(sim_colors.simlist):
    if sim in quand:
        continue
    print("=== %s ==="%sim)
    ol = "%s/%s/OutputLog"%(sim_dir,sim)
    quand[sim] = gaq.get_quantities_and_rms(ol)

#
# Read spectra
#

plotdir = "//home/dccollins/PigPen"
if 'clobber' not in dir():
    clobber=False
if 'spectra_dict' not in dir() or clobber==True:
    spectra_dict={}
    shortprefix='time'
    for axes in ['x','y','z']:
        print('axes',axes)
        spectra_dict[axes]={}
        for i, sim in enumerate(sim_colors.simlist):
            spectra_dict[axes][sim]=slopes(sim) 

            sim_dir = "/data/cb1/Projects/P49_EE_BB/512_frbs/%s"%sim
            frb_name = ""
            longprefix='%s_%s_%s'%(sim,shortprefix,axes)
            spectra_fname = "avg_spectra_%s_%s.h5"%(sim,axes)
            pack = queb3.simulation_package( directory=sim_dir,frames=sim_colors.frames[sim],prefix=longprefix,
                                           frbname=frb_name)
            proj=pack.read_queb(frame=sim_colors.frames[sim][0],ax=axes,bin_style='dx1')
            fitrange=proj.determine_fit_range()  #something better
            if os.path.exists(spectra_fname):
                try:
                    fptr=h5py.File(spectra_fname,'r')
                    for field in fptr:
                        spectra_dict[axes][sim].spectra[field]=fptr[field][()]
                except:
                    raise
                finally:
                    fptr.close()
            spectra_dict[axes][sim].fit(proj.lcent,fitrange)
            spectra_dict[axes][sim].lcent = proj.lcent
            spectra_dict[axes][sim].fit_range = fitrange    


#
# Read restricted set.
#

if 'quan3' not in dir():
    quan3={}
    for ns, sim in enumerate(sim_colors.simlist):
        quan3[sim]={}
        v2avg=[]
        msavg=[]
        maavg=[]
        for frame in sim_colors.framelist[ns]:
            fname = '%s/%s/DD%04d/data%04d.AverageQuantities.h5'%(sim_colors.cloudbreak_base,sim,frame,frame)
            if not os.path.exists(fname):
                print(fname)
                continue
            quan4=h5py.File(fname,'r')
            v2 = np.sqrt(quan4['vx_std'][:]**2+quan4['vy_std'][:]**2+quan4['vz_std'][:]**2)
            bt = np.sqrt(quan4['bx_avg'][:]**2+quan4['by_avg'][:]**2+quan4['bz_avg'][:]**2)
            ma=v2/bt
            v2avg=np.append(v2avg,v2)
            maavg=np.append(maavg,ma)
            print(np.mean(ma), np.mean(maavg))
        quan3[sim]['maavg'] =np.mean(maavg)
        quan3[sim]['msavg'] =np.mean(v2avg)

#
# Spectra plots
#


if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True,figsize=(12,4))
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:
            label="%s %0.1f"%(sim, spectra_dict['x'][sim].slopes[field])
            axlist[nf].plot(proj.lcent,   spectra_dict['x'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)

            label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            #axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    for a in axlist[:3]:
        a.legend(loc=0)
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/alpha_rho_v_H.pdf'%plotdir)


#
# Primitive spectra
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True,figsize=(12,4))
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_d','avg_v','avg_h']):
        for sim in sim_colors.simlist:
            label="%s %0.2f"%(sim, spectra_dict['x'][sim].slopes[field])
            print(label)
            kwargs={'c':sim_colors.color[sim], 'linestyle':sim_colors.linestyle[sim], 
                    'label':label}
            axlist[nf].plot(proj.lcent,   spectra_dict['x'][sim].spectra[field], **kwargs)
            kwargs.pop('label')
            pline(spectra_dict['x'][sim], axlist[nf], field, c='k')
            
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    #for a in axlist[:3]:
    #    a.legend(loc=1)
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[2]]:
        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/spectra_rho_v_H.pdf'%plotdir)

#
# Primitive Slope vs TEB slope
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(3,3, sharex=True,sharey=True,figsize=(12,4))
    #fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nfx,fieldx in enumerate(['avg_d','avg_v','avg_h']):
        for nfy,fieldy in enumerate(['avg_cltt','avg_clee','avg_clbb']):
            for sim in sim_colors.simlist:
                #label="%s %0.2f"%(sim, spectra_dict['x'][sim].slopes[field])
                #print(label)
                #kwargs={'c':sim_colors.color[sim], 'linestyle':sim_colors.linestyle[sim], 
                #        'label':label}
                ax[nfy][nfx].scatter( spectra_dict['y'][sim].slopes[fieldx],
                                    spectra_dict['y'][sim].slopes[fieldy],
                                   c=sim_colors.color[sim], marker=sim_colors.marker[sim])
                #kwargs.pop('label')
                #pline(spectra_dict['x'][sim], axlist[nf], field, c='k')

            
    for a in axlist:
        dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
    ax[0][0].set_ylabel(r'$\alpha_\rho$')
    ax[1][0].set_ylabel(r'$\alpha_v$')
    ax[2][0].set_ylabel(r'$\alpha_H$')
    ax[2][0].set_xlabel(r'$\alpha_{\rm{TT}}$')
    ax[2][1].set_xlabel(r'$\alpha_{\rm{EE}}$')
    ax[2][2].set_xlabel(r'$\alpha_{\rm{BB}}$')


    fig.savefig('%s/slope_TTEEBB_vs_rho_v_H.pdf'%plotdir)
