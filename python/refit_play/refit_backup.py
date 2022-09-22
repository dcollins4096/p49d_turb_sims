
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


#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
#reload(read_stuff)  
spectra_dict = rs.spectra_dict
quan3 = rs.quan3



#
# Spectra plots
#

plotdir=os.environ['HOME']+"/PigPen"
from scipy.stats import ks_2samp
def do_ks_string(test_dist, original_dist):
        KS_output = ks_2samp(test_dist, original_dist)
        crit_stat = 1.36*np.sqrt( (original_dist.size + test_dist.size)/(original_dist.size*test_dist.size))
        ks_string = r'$D=%0.2f (D_c=%0.2f, p=%f)$'%(KS_output.statistic,crit_stat, KS_output.pvalue)
        return ks_string
def do_ks(test_dist, original_dist):
        KS_output = ks_2samp(test_dist, original_dist)
        crit_stat = 1.36*np.sqrt( (original_dist.size + test_dist.size)/(original_dist.size*test_dist.size))
        #ks_string = r'$D=%0.2f (D_c=%0.2f, p=%f)$'%(KS_output.statistic,crit_stat, KS_output.pvalue)
        return KS_output.statistic,crit_stat, KS_output.pvalue

if 1:
    plt.close('all')
    #fig,ax = plt.subplots(2,3, sharex=True,sharey=True,figsize=(12,4))
    fig,ax = plt.subplots(1,1, sharex=True,sharey=True,figsize=(12,4))
    fig.subplots_adjust(wspace=0, hspace=0)
    #axlist=ax.flatten()
    axlist=[ax]

    counter=0
    slope_old=[]
    slope_new=[]    
    f3,ax3=plt.subplots(1,1)
    #ax2=ax2.flatten()

    #for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb'])
    chiext=dt.extents()
    slopext=dt.extents()
    for nf,field in enumerate(['avg_clbb']):
        #for aaa in ax2:
        #    aaa.clear()
        ax3.clear()
        for ns,sim in enumerate(sim_colors.simlist):
            f2,ax2=plt.subplots(2,2)
            ax2=ax2.flatten()
            counter += 1
            print(counter)
            #if counter not in [2]:
            #    continue
            #if sim != '1_2':
            #    continue

            proj=rs.proj_dict['x'][sim]
            label="%s %0.1f"%(sim, spectra_dict['x'][sim].slopes[field])

            x0_list = proj.lcent[2:10:2]
            x1_list = proj.lcent[11:41:4]
            max_npoints=np.argmax(x1_list)-np.argmin(x0_list)+1
            rm = dt.rainbow_map(max_npoints)
            rm0 = dt.rainbow_map(len(x0_list))
            rm = dt.rainbow_map(len(x1_list))

            fig4,ax4=plt.subplots(1,1)
            ntot=0
            spectra_dict[LOS][sim].spext=dt.extents()
            spectra_dict[LOS][sim].speyt=dt.extents()
            xlim = {'1_half':[3.07e-02, 4.85e-01], '1_2':[3.07e-02, 4.85e-01]}.get(sim,None)

            ylim = {'1_half':[1.77e-05, 5.93e-01], '1_2':[4.73e-03, 1.03e+02]}.get(sim,None)
            slopelim=  [-4.75e+00, -1.81e+00]
            chilim = [1.24e-04, 2.05e+03]

            spectra_dict[LOS][sim].slope_array= np.zeros([x0_list.size, x1_list.size]) +10
            spectra_dict[LOS][sim].chi_array= np.zeros([x0_list.size, x1_list.size]) +10
            for n0,x0 in enumerate(x0_list):
                for n1,x1 in enumerate(x1_list):

                    ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                    xvals = proj.lcent[ok]
                        
                    #axlist[nf].plot(proj.lcent[ok],   spectra_dict['x'][sim].spectra[field][ok], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)

                    spec=spectra_dict['x'][sim].spectra[field][ok]
                    slope,amp,res=plfit(xvals,spec,[x0,x1])

                    Line = 10**( slope*np.log10(xvals)+np.log10(amp))
                    #axlist[nf].plot( xvals, Line,c='k',label='F1')
                    spectra_dict[LOS][sim].slope_array[n0,n1] = slope
                    if slope>0:
                        continue    

                    #KS_stat, KS_crit, KS_p=do_ks( np.log10(Line), np.log10(spec))
                    ##ax2.scatter( slope, ( (Line-spec)**2/spec).sum(), c=[rm( ok.sum())])
                    #chi_2 = ((Line-spec)**2/spec).sum()/( ok.sum()-2)
                    chi_2 = ((np.log10(Line)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)

                    spectra_dict[LOS][sim].chi_array[n0,n1] = chi_2
                    spectra_dict[LOS][sim].slopext(nar([slope]))
                    spectra_dict[LOS][sim].chiext(nar([chi_2]))
                    #continue
                    ax2[0].scatter( slope, chi_2, c=[rm( ok.sum())])
                    #ax2[0].scatter( slope, chi_2, c=[rm(n0)])
                    ax2[2].scatter( ok.sum(), chi_2,c='k')
                    #ax2[2].scatter( ok.sum(), chi_2/ok.sum(), c='r')
                    #ax2[3].scatter( ok.sum(), slope, c=[rm(n1)])
                    ax2[3].scatter( ok.sum(), slope, c=[rm0(n0)])

                    ax4.clear()
                    ax4.plot( xvals, Line, c='g')
                    ax4.plot( xvals, spec, c='k')
                    spext(xvals);speyt(Line); speyt(spec)
                    dt.axbonk(ax4,xlim=xlim, ylim=ylim,xscale='log',yscale='log')
                    fig4.savefig('%s/OOF_%s_%d.png'%(plotdir,sim,ntot))
                    ntot+=1 


            if 0:
                n = spec.spectra[field].size
                alphahat = 1 + ok.sum()/( np.log(spec.spectra[field][ok]/fitrange[0]).sum())
                print("SLOOOO", spec.slopes[field], alphahat)
                dk = proj.lbins[1:]-proj.lbins[:-1]
                A = np.sum( spec.spectra[field][ok]*dk[ok])*(1-alphahat)*fitrange[0]**alphahat/(fitrange[1]**(1-alphahat)-fitrange[0]**(1-alphahat))

                Line2 = 10**( alphahat*np.log10(proj.lcent))
                Line2 *= spec.spectra[field][ok][5]/Line2[ok][5]
                axlist[nf].plot( proj.lcent[ok], Line2[ok],label='F2',c='r')
            #bins = np.mgrid[ slope_array.min():slope_array.max():16j]
            bins = np.mgrid[ -4.5:0:16j] +np.random.random()*0.01   
            bins = np.mgrid[slope_array.min(): slope_array.min()+4.5:0.3]
            print(bins)
            #if slope_array.min()<bins.min():
            #    raise
            #if slope_array.max()>bins.max():
            #    raise

            #ax2[1].hist( slope_array[slope_array < 0].flatten())
            ax3.clear()
            ax3.hist( slope_array[slope_array < 0].flatten(), histtype='step', linestyle=sim_colors.line_list[ns],color=sim_colors.color_list[ns],bins=bins)
            dt.axbonk(ax3,xlabel=field,ylabel='N',xlim=[-4.5,0],ylim=[0,200])
            f3.savefig('%s/refit_SlopeHist_%s_%s.png'%(plotdir,field,sim))

            ax2[1].hist( slope_array[slope_array < 0].flatten(), histtype='step', linestyle=sim_colors.line_list[ns],color=sim_colors.color_list[ns],bins=bins)
            dt.axbonk(ax2[1],xlabel=field,ylabel='N',xlim=[-4.5,0],ylim=[0,200])
            dt.axbonk(ax2[0],xlabel='slope',ylabel='Chi',yscale='log', ylim=chilim,xlim=slopelim)
            dt.axbonk(ax2[2],xlabel='Ndof',ylabel='Chi',yscale='log', ylim=chilim)
            dt.axbonk(ax2[3],xlabel='Ndof',ylabel='Slope', ylim=slopelim)


            #proj=rs.proj_dict['y'][sim]
            #label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            #axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            #axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            f2.savefig('%s/refit_SlopeThings_%s_%s.png'%(plotdir,field,sim))
            plt.close(f2)
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    for a in axlist[:3]:
        a.legend(loc=0)
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
#    for a in [axlist[0],axlist[3]]:
#        a.set_ylabel(r'$P(k)$')

    fig.savefig('%s/refit_eb.pdf'%plotdir)


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
