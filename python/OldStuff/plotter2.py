
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
if 0:
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

#
# Prim amps nominal, old code 
#

if 0:
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


#
# Measured vs. nominal Ms Ma
#


if 0:
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
            

