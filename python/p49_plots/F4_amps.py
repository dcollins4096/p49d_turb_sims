
from GL import *

import simulation
import simulation_info.all_sims

def plot_amps_slopes(sim_list, prim_or_teb='prim',amps_or_slopes='amps', axis='y', fit_herd=None):

    do_prim=False;do_TEB=False;do_amp=False;do_slope=False;do_log=False
    if prim_or_teb == 'prim':
        product_list = ['density','velocity','Htotal']
        axis_label=""
        do_prim=True
        S1,S2,S3=r'\rho','v','H'
    else:
        product_list = ['ClTT'+axis,'ClEE'+axis,'ClBB'+axis]
        axis_label="_%s"%axis
        do_TEB=True
        S1,S2,S3='TT','EE','BB'
    if amps_or_slopes=='amps':
        do_amp = True
        Aalpha='A'
        do_log=True
    else:
        do_slope=True
        Aalpha='\\alpha'

    fig,axes = plt.subplots(3,2, figsize=(5.5,5.5))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=axes.flatten()

    ext = dt.extents()
    for ns, sim in enumerate(sim_list):
        this_sim = simulation.corral[sim]
        this_sim.load()
        #simulation.set_colors(this_sim)
        simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)
        ms = this_sim.Ms_mean
        ma = this_sim.Ma_mean
        for nf, field in enumerate(product_list):
            if amps_or_slopes == 'amps':
                y = this_sim.ampsA[field]
            elif amps_or_slopes == 'slopes':
                y = this_sim.slopesA[field]
            axes[nf][0].scatter(ms, y, c=[this_sim.color],marker=this_sim.marker, s=this_sim.marker_size*20)
            axes[nf][1].scatter(ma, y, c=[this_sim.color],marker=this_sim.marker, s=this_sim.marker_size*20)
            if do_log and y<0:
                pdb.set_trace()
            ext(y)

    if do_log:
        for aaa in axlist:
            aaa.set_yscale('log')
    axes[0][0].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S1))
    axes[1][0].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S2))
    axes[2][0].set_ylabel(r'$%s_{\rm{%s}}$'%(Aalpha,S3))
    axes[2][0].set_xlabel(r'$M_{\rm{s}}$')
    axes[2][1].set_xlabel(r'$M_{\rm{A}}$')

    if (do_prim and do_amp):
        for aaa in axlist:
            #aaa.set_ylim([4e-6,4e-1])
            aaa.set_ylim([5e-7,2e-2])
            #print('word')
    if (do_TEB and do_amp):
        for aaa in axlist:
            aaa.set_ylim([2e-7,2e-1])
    if do_prim and do_slope:
        for aaa in axlist:
            pass
            #aaa.set_ylim([-1.9,-0.9])
    if do_TEB and do_slope:
        for aaa in axlist:
            pass
            #aaa.set_ylim([-5,-1])

    if do_log:
        for aaa in axlist:
            y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
            aaa.yaxis.set_minor_locator(y_minor)
            aaa.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        for aaa in axlist:
            aaa.minorticks_on()

    for n in range(3):
        axes[n][1].yaxis.tick_right()

    fig.subplots_adjust(top=0.92)
    #fig.tight_layout()
    if amps_or_slopes=='amps':
        fig.suptitle('Amplitudes', y=0.96)
    else:
        fig.suptitle('Slopes', y=0.96)


    fit_text=""
    if fit_herd is not None:
        fit_text="with_fit"
        for nf, field in enumerate(product_list):
            if amps_or_slopes=='amps':
                this_prod = field+"_a"
            else:
                this_prod = field+"_s"
            fff = fit_herd[this_prod]
            params=fff.Params
            if len(params)==3:
                a,b,c=params
            else:
                a,b,c,d=params
            thax=axes[nf][0]
            msarr=nar([0.5,1,2,3,4,5,6])
            #this_ma = fff.Ma.mean()
            #print(fff.Ma, this_ma)
            #this_ma=0.5
            
            for nma,ma in enumerate([0.5,1,2]):
                if len(params)==4:
                    line = a + msarr*b + ma*c+ma*msarr*d
                else:
                    line = a + msarr*b + ma*c
                color='rgb'[nma]
                thax.plot(msarr,line, c=color)

        


    axislabel=''
    if prim_or_teb=='teb':
        axislabel='_%s'%axis
    outname = "%s/AS_%s%s_%s%s.pdf"%(dl.plotdir,prim_or_teb,axislabel,amps_or_slopes,fit_text)
    fig.savefig(outname)
    print(outname)

