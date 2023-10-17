
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

#
# Read restricted set.
#

if 0:
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

    fig.savefig('%s/TT_EE_BB.pdf'%plotdir)

#
# TT EE BB vs Mach
#

if 0:
    for xyz in 'xyz':
        plt.close('all')
        fig,ax = plt.subplots(3,2)#,sharey=True)
        fig.subplots_adjust(wspace=0, hspace=0)
        for sim in sim_colors.simlist:
            #nominal or measured
            if 1:
                this_Ms = quand[sim]['Ms'].mean()
                this_Ma = quand[sim]['Ma'].mean()
                measured_or_nominal = 'measured'
            else:
                this_Ms = sim_colors.Ms[sim]
                this_Ma = sim_colors.Ma[sim]
                measured_or_nominal = 'nominal'

            if 1:
                amps_or_slopes = spectra_dict[xyz][sim].amps
                aos = 'amps'
                Aalpha = 'A'
                yscale='log'
            else:
                amps_or_slopes = spectra_dict[xyz][sim].slopes
                aos = 'slopes'
                Aalpha = "\\alpha"
                yscale='linear'
            outname = '%s/TTEEBB_%s_mach_%s_%s.pdf'%(plotdir, xyz, aos, measured_or_nominal)

            label="%s %0.1f"%(sim, spectra_dict[xyz][sim].slopes['avg_cltt'])
            ax[0][0].scatter(this_Ma, amps_or_slopes['avg_cltt'], c=sim_colors.color[sim], marker=sim_colors.marker[sim],label=label)
            ax[1][0].scatter(this_Ma, amps_or_slopes['avg_clee'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[2][0].scatter(this_Ma, amps_or_slopes['avg_clbb'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            
            ax[0][1].scatter(this_Ms, amps_or_slopes['avg_cltt'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[1][1].scatter(this_Ms, amps_or_slopes['avg_clee'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[2][1].scatter(this_Ms, amps_or_slopes['avg_clbb'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])

            ax[0][0].set_ylabel(r'$%s_{TT}$'%Aalpha)
            ax[1][0].set_ylabel(r'$%s_{EE}$'%Aalpha)
            ax[2][0].set_ylabel(r'$%s_{BB}$'%Aalpha)
            ax[2][0].set_xlabel(r'$M_{\rm{A}}$')
            ax[2][1].set_xlabel(r'$M_{\rm{s}}$')
            for a in [ax[0][1],ax[1][1],ax[2][1]]:
                a.yaxis.tick_right()
        ax[0][0].legend(loc=0)
        for a in ax.flatten():
            a.set_yscale(yscale)
        fig.savefig(outname)
        print(outname)


#
# Primitives vs Mach
#
if 1:
    plt.close('all')
    fig,ax = plt.subplots(3,3)#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for sim in sim_colors.simlist:
        #nominal or measured
        if 1:
            this_Ms = quan3[sim]['msavg'].mean()
            this_Ma = quan3[sim]['maavg'].mean()
            measured_or_nominal = 'measured3'
        elif 0:
            this_Ms = quand[sim]['Ms'].mean()
            this_Ma = quand[sim]['Ma'].mean()
            measured_or_nominal = 'measured'
        else:
            this_Ms = sim_colors.Ms[sim]
            this_Ma = sim_colors.Ma[sim]
            measured_or_nominal = 'nominal'

        if 1:
            amps_or_slopes = spectra_dict['y'][sim].amps
            aos = 'amps'
            Aalpha = 'A'
        else:
            amps_or_slopes = spectra_dict['y'][sim].slopes
            aos = 'slopes'
            Aalpha = "\\alpha"

        MeanB = (quand[sim]['Btot']**2/(8*np.pi)).mean()
        outname = '%s/prim_mach_%s_%s.pdf'%(plotdir, aos, measured_or_nominal)

        ax[0][0].scatter(this_Ma, amps_or_slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][0].scatter(this_Ma, amps_or_slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][0].scatter(this_Ma, amps_or_slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        
        ax[0][1].scatter(this_Ms, amps_or_slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][1].scatter(this_Ms, amps_or_slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][1].scatter(this_Ms, amps_or_slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])

        ax[0][2].scatter(MeanB, amps_or_slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][2].scatter(MeanB, amps_or_slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][2].scatter(MeanB, amps_or_slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])

        ax[0][0].set_ylabel(r'$%s_v$'%Aalpha)
        ax[1][0].set_ylabel(r'$%s_\rho$'%Aalpha)
        ax[2][0].set_ylabel(r'$%s_H$'%Aalpha)
        ax[2][0].set_xlabel(r'$M_{\rm{A}}$')
        ax[2][1].set_xlabel(r'$M_{\rm{s}}$')
        for a in [ax[0][1],ax[1][1],ax[2][1]]:
            a.yaxis.tick_right()
    fig.savefig(outname)
    print(outname)


#
# reb, ttb, rte spectra
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_rte', 'avg_rtb', 'avg_reb']):
        for sim in sim_colors.simlist:
            qu = field[-2:].upper()
            label = r'$r_{%s} \parallel$'%qu
            axlist[nf  ].text(0.1, -0.8, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf+3].text(0.1, -0.8, label)
            axlist[nf  ].plot(proj.lcent, spectra_dict['x'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            #axlist[nf+6].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        a.set_yscale('symlog',linthresh=0.09)
        #a.set_yscale('linear')
        a.set_ylim([-1,1])
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$r_{XY}$')
        a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    fig.savefig('%s/r_XY.pdf'%plotdir)

#
# reb, ttb, rte PDFs
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()
    bins = np.linspace(-1,1,64)
    linear = 0.1
    b1 = np.linspace(-linear, linear, 8)
    b2 = np.logspace(np.log10(linear), 0, 8)
    bins = np.concatenate([-b2[::-1], b1, b2])
    bins = np.unique(bins)
    for sim in sim_colors.simlist:
        #for a in ax.flatten():
        #    a.clear()
        for nf,field in enumerate(['avg_rte', 'avg_rtb', 'avg_reb']):
            qu = field[-2:].upper()
            label = r'$r_{%s} \parallel$'%qu
            axlist[nf  ].text(-0.5, 100, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf+3].text(-0.5, 100, label)
            h1,b1=np.histogram(spectra_dict['x'][sim].spectra[field],bins=bins)
            h2,b2=np.histogram(spectra_dict['y'][sim].spectra[field],bins=bins)
            b1c = 0.5*(b1[:-1]+b1[1:])
            b2c = 0.5*(b2[:-1]+b2[1:])
            b1w = (b1[1:]-b1[:-1])
            b2w = (b2[1:]-b2[:-1])

            axlist[nf  ].plot(b1c, h1, color=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            axlist[nf+3].plot(b2c, h2, color=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
        #plt.savefig("%s/r_XY_pdf_%s.pdf"%(plotdir,sim))
    for a in axlist:
        dt.axbonk(a,yscale='log',xlabel=None,ylabel=None)
        a.set_xscale('symlog',linthresh=0.151515151515151515151515151515)
#        a.set_ylim([-1,1])
    for a in axlist[3:]:
        a.set_xlabel(r'$r_{XY}$')
#   for a in [axlist[0],axlist[3]]:
#       a.set_ylabel(r'$P(r_{XY})$')
#       a.set_yticks([-1,-0.1,-0.01,0.01,0.1,1])

    fig.savefig('%s/r_XY_pdf.pdf'%plotdir)

#
# reb, ttb, rte mean and variance
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()
    for sim in sim_colors.simlist:
        #for a in ax.flatten():
        #    a.clear()
        for nf,field in enumerate(['avg_rte', 'avg_rtb', 'avg_reb']):
            qu = field[-2:].upper()
            label = r'$r_{%s} \parallel$'%qu
            axlist[nf  ].text(0.5e-3, 0.2, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf+3].text(0.5e-3, 0.2, label)

            mean1 = spectra_dict['x'][sim].spectra[field].mean()
            std1  = spectra_dict['x'][sim].spectra[field].std()
            mean2 = spectra_dict['y'][sim].spectra[field].mean()
            std2  = spectra_dict['y'][sim].spectra[field].std()

            axlist[nf  ].scatter(np.abs(mean1),std1, color=sim_colors.color[sim], marker=sim_colors.marker[sim])
            axlist[nf+3].scatter(np.abs(mean2),std2, color=sim_colors.color[sim], marker=sim_colors.marker[sim])

    for a in axlist:
        #dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    for a in axlist[3:]:
        a.set_xlabel(r'$\langle r_{XY} \rangle$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$\sigma_{XY}$')

    fig.savefig('%s/r_XY_mean_std.pdf'%plotdir)



if 0:
    for xyz in 'xyz':
        plt.close('all')
        fig,ax = plt.subplots(3,3,sharey=True, sharex=True)
        axlist=ax.flatten()
        fig.subplots_adjust(wspace=0, hspace=0)
        for sim in sim_colors.simlist:
            
            ax[0][0].set_ylabel(r'$\alpha_v$')
            ax[1][0].set_ylabel(r'$\alpha_\rho$')
            ax[2][0].set_ylabel(r'$\alpha_H$')
            ax[2][0].set_xlabel(r'$C_\ell^{TT}$')
            ax[2][1].set_xlabel(r'$C_\ell^{EE}$')
            ax[2][2].set_xlabel(r'$C_\ell^{BB}$')

            ax[0][0].scatter(spectra_dict[xyz][sim].amps['avg_cltt'], spectra_dict[xyz][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[1][0].scatter(spectra_dict[xyz][sim].amps['avg_cltt'], spectra_dict[xyz][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[2][0].scatter(spectra_dict[xyz][sim].amps['avg_cltt'], spectra_dict[xyz][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])

            ax[0][1].scatter(spectra_dict[xyz][sim].amps['avg_clee'], spectra_dict[xyz][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[1][1].scatter(spectra_dict[xyz][sim].amps['avg_clee'], spectra_dict[xyz][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[2][1].scatter(spectra_dict[xyz][sim].amps['avg_clee'], spectra_dict[xyz][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])

            ax[0][2].scatter(spectra_dict[xyz][sim].amps['avg_clbb'], spectra_dict[xyz][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[1][2].scatter(spectra_dict[xyz][sim].amps['avg_clbb'], spectra_dict[xyz][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            ax[2][2].scatter(spectra_dict[xyz][sim].amps['avg_clbb'], spectra_dict[xyz][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            for a in [ax[0][2],ax[1][2],ax[2][2]]:
                a.yaxis.tick_right()
            for val in [0,1,2]:
                ax[val][1].set_ylim( ax[val][0].get_ylim())
                ax[val][2].set_ylim( ax[val][0].get_ylim())
        for a in axlist:
            a.set_xscale('log')
        fig.savefig('%s/ee_bb_prim_amps_%s.pdf'%(plotdir,xyz))

if 0:
    plt.close('all')
    fig,ax = plt.subplots(3,3,sharey=True, sharex=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for sim in sim_colors.simlist:
        ax[0][0].scatter(spectra_dict['y'][sim].slopes['avg_cltt'], spectra_dict['y'][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][0].scatter(spectra_dict['y'][sim].slopes['avg_cltt'], spectra_dict['y'][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][0].scatter(spectra_dict['y'][sim].slopes['avg_cltt'], spectra_dict['y'][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])

        ax[0][1].scatter(spectra_dict['y'][sim].slopes['avg_clee'], spectra_dict['y'][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][1].scatter(spectra_dict['y'][sim].slopes['avg_clee'], spectra_dict['y'][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][1].scatter(spectra_dict['y'][sim].slopes['avg_clee'], spectra_dict['y'][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        print(sim_colors.Ma[sim], spectra_dict['y'][sim].slopes['avg_d'])
        
        ax[0][0].set_ylabel(r'$\alpha_v$')
        ax[1][0].set_ylabel(r'$\alpha_\rho$')
        ax[2][0].set_ylabel(r'$\alpha_H$')
        ax[2][0].set_xlabel(r'$C_\ell^{TT}$')
        ax[2][1].set_xlabel(r'$C_\ell^{EE}$')
        ax[2][2].set_xlabel(r'$C_\ell^{BB}$')
        ax[0][2].scatter(spectra_dict['y'][sim].slopes['avg_clbb'], spectra_dict['y'][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][2].scatter(spectra_dict['y'][sim].slopes['avg_clbb'], spectra_dict['y'][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][2].scatter(spectra_dict['y'][sim].slopes['avg_clbb'], spectra_dict['y'][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        for a in [ax[0][2],ax[1][2],ax[2][2]]:
            a.yaxis.tick_right()
        for val in [0,1,2]:
            ax[val][1].set_ylim( ax[val][0].get_ylim())
            ax[val][2].set_ylim( ax[val][0].get_ylim())
    fig.savefig('%s/ee_bb_prim_spectra.pdf'%plotdir)


if 1:
    plt.close('all')
    fig,ax = plt.subplots(3,2)#,sharey=True)

    fig.subplots_adjust(wspace=0, hspace=0)
    for sim in sim_colors.simlist:
        ax[0][0].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].amps['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][0].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].amps['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][0].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].amps['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        
        ax[0][0].set_ylabel(r'$\alpha_v$')
        ax[1][0].set_ylabel(r'$\alpha_\rho$')
        ax[2][0].set_ylabel(r'$\alpha_H$')
        ax[2][0].set_xlabel(r'$M_{\rm{A}}$')
        ax[2][1].set_xlabel(r'$M_{\rm{s}}$')
        ax[0][1].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].amps['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][1].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].amps['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][1].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].amps['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        for a in [ax[0][1],ax[1][1],ax[2][1]]:
            a.yaxis.tick_right()
    for a in ax.flatten():
        a.set_yscale('log')
    fig.savefig('%s/prim_amps.pdf'%plotdir)



if 1:
    plt.close('all')
    fig,ax = plt.subplots(2,2,sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for i,sim in enumerate(sim_colors.simlist):
        this_Ms = quand[sim]['Ms'].mean() 
        this_Ma = quand[sim]['Ma'].mean() 
        this_Ma2,this_Ms2=sim_colors.Ma[sim], sim_colors.Ms[sim]
        ax[0][1].scatter(this_Ma, this_Ms, c=sim_colors.color[sim], marker=sim_colors.marker[sim],label=sim)
        ax[0][1].plot( [this_Ma2, this_Ma], [this_Ms2,  this_Ms], c=sim_colors.color[sim])
        ax[0][0].scatter(this_Ma2, this_Ms2, c=sim_colors.color[sim], marker=sim_colors.marker[sim],label=sim)
        #ax.plot([this_Ma,this_Ma2],[this_Ms,this_Ms2])
        ax[1][0].scatter( alfmachavg[i],sonicmachavg[i])
        ax[1][1].scatter( quan3[sim]['maavg'], quan3[sim]['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][1].plot( [this_Ma2, quan3[sim]['maavg']], [this_Ms2,  quan3[sim]['msavg']], c=sim_colors.color[sim])
    dt.axbonk(ax[0][0],xlabel=r'$M_{\rm{A}}$',ylabel=r'$M_{\rm{S}}$')#,xlim=[0.45,2.1],ylim=[0.45,3.2])
    fig.savefig('%s/point_legend_measured.pdf'%plotdir)
            

#
# Primitives Lines
#
AAA=defaultdict(list)
MMMa=defaultdict(list)
MMMs=defaultdict(list)
matrix_d = np.zeros([4,3]) #ms,ma
matrix_v = np.zeros([4,3]) #ms,ma
matrix_h = np.zeros([4,3]) #ms,ma
the_x = np.zeros_like(matrix)
the_y = np.zeros_like(matrix)
if 1:
    plt.close('all')
    fig,ax = plt.subplots(3,3, figsize=(8,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axis = 'y'
    simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
    sublist=nar([["half_half","half_1","half_2"],["1_half","1_1","1_2"],["2_half","2_1","2_2"],["3_half","3_1","3_2"]]).transpose()

    for nsub,sub  in enumerate(sublist):
        #nominal or measured
        if 1:
            this_Ms =[ quan3[sim]['msavg'].mean() for sim in sub]
            this_Ma =[ quan3[sim]['maavg'].mean() for sim in sub]
            measured_or_nominal = 'measured3'
        elif 0:
            this_Ms =[ quand[sim]['Ms'].mean() for sim in sub]
            this_Ma =[ quand[sim]['Ma'].mean() for sim in sub]
            measured_or_nominal = 'measured'
        else:
            this_Ms = [sim_colors.Ms[sim]]*len(sub)
            this_Ma = [sim_colors.Ma[sim]]*len(sub)
            measured_or_nominal = 'nominal'


        if 1:
            avg_v  = [spectra_dict[axis][sim].amps['avg_v'] for sim in sub]
            avg_d  = [spectra_dict[axis][sim].amps['avg_d'] for sim in sub]
            avg_h  = [spectra_dict[axis][sim].amps['avg_h'] for sim in sub]
            aos = 'amps'
            Aalpha = 'A'
        else:
            avg_v  = [spectra_dict[axis][sim].slopes['avg_v'] for sim in sim_colors.simlist]
            avg_d  = [spectra_dict[axis][sim].slopes['avg_d'] for sim in sim_colors.simlist]
            avg_h  = [spectra_dict[axis][sim].slopes['avg_h'] for sim in sim_colors.simlist]
            aos = 'slopes'
            Aalpha = "\\alpha"


        MeanB = nar([(quand[sim]['Btot']**2/(8*np.pi)).mean() for sim in sub])
        outname = '%s/prim_mach_%s_%s.pdf'%(plotdir, aos, measured_or_nominal)
        colors =[sim_colors.color[sim] for sim in sub]
        print(colors)
        markers=[sim_colors.marker[sim] for sim in sub]

        #Its dumb as hell that marker isn't a list
        #and I have to make this loop
        for i in range(len(markers)):
            ax[0][0].scatter(this_Ma[i], avg_v[i], c=colors[i],marker=markers[i])
            ax[1][0].scatter(this_Ma[i], avg_d[i], c=colors[i],marker=markers[i])
            ax[2][0].scatter(this_Ma[i], avg_h[i], c=colors[i],marker=markers[i])

            ax[0][1].scatter(this_Ms[i], avg_v[i], c=colors[i],marker=markers[i])
            ax[1][1].scatter(this_Ms[i], avg_d[i], c=colors[i],marker=markers[i])
            ax[2][1].scatter(this_Ms[i], avg_h[i], c=colors[i],marker=markers[i])

            ax[0][2].scatter(  MeanB[i], avg_v[i], c=colors[i],marker=markers[i])
            ax[1][2].scatter(  MeanB[i], avg_d[i], c=colors[i],marker=markers[i])
            ax[2][2].scatter(  MeanB[i], avg_h[i], c=colors[i],marker=markers[i])

        fits_v = np.polyfit(MeanB, avg_v, 1)
        fits_d = np.polyfit(MeanB, avg_d, 1)
        fits_h = np.polyfit(MeanB, avg_h, 1)
        #ax[0][2].plot(MeanB, fits_v[0]*MeanB+fits_v[1], c=colors[0])
        #ax[1][2].plot(MeanB, fits_d[0]*MeanB+fits_d[1], c=colors[0])
        #ax[2][2].plot(MeanB, fits_h[0]*MeanB+fits_h[1], c=colors[0])

        matrix_d[:,nsub]=avg_d
        matrix_v[:,nsub]=avg_v
        matrix_h[:,nsub]=avg_h
        the_x[:,nsub]=MeanB
        the_y[:,nsub]=this_Ms


        ax[0][0].set_ylabel(r'$%s_v$'%Aalpha)
        ax[1][0].set_ylabel(r'$%s_\rho$'%Aalpha)
        ax[2][0].set_ylabel(r'$%s_H$'%Aalpha)
        ax[2][0].set_xlabel(r'$M_{\rm{A}}$')
        ax[2][1].set_xlabel(r'$M_{\rm{s}}$')

    def do_stuff(xdata, ydata, p0):
        def fitFunc(x, a, b, c, d):
            return a+ b*x[0] + c*x[1] + d*x[0]*x[1]
        fitParams, fitCovariances = curve_fit(fitFunc, xdata, ydata, p0)
        print(' fit coefficients:\n', fitParams)
        return fitParams, fitCovariances
    def do_morestuff(aax, myx, myy, myval,c='k'):
        x=np.row_stack([myx.flatten(),myy.flatten()])
        y=myval.flatten()
        p0 = [0,1,1,0]
        fitParams, fitCovariances= do_stuff(x,y, p0)
        f = fitParams
        amps = f[0] + f[1]*myx + f[2]*myy+f[3]*myx*myy
        for b, sim in enumerate(sim_colors.simlist):
            aax.plot( myx[:,b], amps[:,b],c=c[b])
            aax.scatter( myx[:,b], amps[:,b], c=sim_colors.color[sim])# myx[:,nsub], amps,c='k')
            print('uuu',c)
        #ax[1][2].scatter(the_x, matrix_d,c='k')
    do_morestuff( ax[0][2], the_x, the_y, matrix_v,c=colors)
    do_morestuff( ax[1][2], the_x, the_y, matrix_d,c=colors)
    do_morestuff( ax[2][2], the_x, the_y, matrix_h,c=colors)

    for a in [ax[0][1],ax[1][1],ax[2][1]]:
        a.set_yticks([])
    for a in [ax[0][-1],ax[1][-1],ax[2][-1]]:
        a.yaxis.tick_right()

    fig.savefig(outname)
    print(outname)

if 0:
    nominal_Ma =nar([sim_colors.Ma[sim] for sim in sub])
    nominal_Ms =nar([sim_colors.Ms[sim] for sim in sub])
    ax[0][3].scatter(nominal_Ma.mean(), fits_v[0], c=colors[0])
    AAA['v'].append(fits_v[0])
    MMMa['v'].append(nominal_Ma.mean())
    MMMa['v']=nar(MMMa['v'])
    AAA['v']=nar(AAA['v'])
    pfit = np.polyfit(MMMa['v'],AAA['v'],1)
    ax[0][3].plot(MMMa['v'],pfit[0]*MMMa['v']+pfit[1],c='k')

