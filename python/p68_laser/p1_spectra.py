from GL import *
import simulation
import simulation_info.all_sims as all_sims

import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
reload(bt)

import bucket
def brunt_spectra(simlist, ftool=None):
    plt.close('all')

    ncol = int(np.ceil(len(simlist)/3))
    fig,axes = plt.subplots(ncol,3,figsize=(12,12))
    if len(axes)==1:
        axlist=[axes]
    else:
        axlist=axes.flatten()

    #fig.subplots_adjust(wspace=0, hspace=0)

    f2,axes2=plt.subplots(2,2)
    a2=axes2[0][0];a2b=axes2[0][1];a2c=axes2[1][0];a2d=axes2[1][1]
    for i,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]
        this_sim.load()
        simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)                                                                                                                 #ax.scatter( raq.quan3[sim]['maavg'], raq.quan3[sim]['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim],s=60)
        frame=this_sim.ann_frames[-1]


        
        ftool=bucket.things.get(sim,None)
        if ftool is None:
            rho_full, rho = bt.get_cubes(sim,frame)

            ftool = bt.fft_tool(rho)
            ftool.do3()
            ftool.do2(projax=0)
        bucket.things[sim]=ftool

        ax=axlist[i]


        fit_range =this_sim.get_fitrange(xvals)
        mask = (xvals > fit_range[0])*(xvals < fit_range[1])
        ax.plot(ftool.k2d,          ftool.power_1d2.real,c=[0.5]*3, label='P2d')
        ax.plot(ftool.k2d[mask],          ftool.power_1d2.real[mask],c='r', label='P2d')
        ax.plot(ftool.k3d[mask],          ftool.power_1d3.real[mask],c='g', label='P3d')
        ax.plot(ftool.k2d[mask],ftool.k2d[mask]*ftool.power_1d2.real[mask],c='b', label = 'k P2d')
        ax.set(xscale='log',yscale='log')
        #pdb.set_trace()


        if 0:
            sigma_rho = (ftool.rho**2).sum().real
            sigma_fft = (ftool.power_1d3).sum().real
            #sigma_fft = (ftool.power).sum().real/ftool.power.size**2
            #print("rho %0.2e fft %0.2e ratio %0.2e"%(sigma_rho,sigma_fft,(sigma_fft)/sigma_rho))

            sigma_col = ((ftool.rho2)**2).sum().real
            sigma_2d_fft = ftool.power_1d2.sum().real
            #print("col %0.2e fft %0.2e ratio %0.2e"%(sigma_col, sigma_2d_fft, sigma_col/sigma_2d_fft))
            Rinv =( ( ftool.k2d*ftool.power_1d2)[1:].sum()/(ftool.power_1d2)[1:].sum()).real
            Rinv_actual = ftool.power_1d3.sum()/ftool.power_1d2.sum()
            sigma_Brunt = sigma_col.real*Rinv
            sigma_Brunt_actual = sigma_col*Rinv_actual
            #print( sigma_Brunt/sigma_rho)
            #print( sigma_Brunt_actual/sigma_rho)
            R1 =sigma_Brunt/sigma_rho
            R2 = sigma_Brunt_actual/sigma_rho
            #ax.set(title='%s %0.1e %0.1e'%(sim,R2,R1))
            ax.text(0.1,0.1,"%s"%sim,transform=ax.transAxes)
            ax.text(0.25,0.1,"%0.1f %0.1f"%(R2,R1),transform=ax.transAxes)

        if 1:
            #works pretty good
            Nz = ftool.rho.size
            mean_rho=(ftool.rho).sum()/Nz
            sigma_rho = ((ftool.rho-mean_rho)**2).sum().real/Nz
            sigma_fft = (ftool.power_1d3[1:]).sum().real/Nz
            sigma_fake_fft = ( ftool.k2d*ftool.power_1d2)[1:].sum()/Nz
            #sigma_fft = (ftool.power).sum().real/ftool.power.size**2
            #print("rho   %0.2e fft %0.2e ratio %0.2e"%(sigma_rho,sigma_fft,(sigma_fft)/sigma_rho))
            #print("fake  %0.2e"%sigma_fake_fft)

            N2d = ftool.rho2.size
            mean_column = ftool.rho2.sum()/N2d
            sigma_col = ((ftool.rho2-mean_column)**2).sum().real/N2d
            sigma_2d_fft = ftool.power_1d2[1:].sum().real/N2d
            #print("col %0.2e fft %0.2e ratio %0.2e"%(sigma_col, sigma_2d_fft, sigma_col/sigma_2d_fft))
            Rinv = (( ftool.k2d*ftool.power_1d2)[1:].sum()/Nz)/((ftool.power_1d2)[1:].sum()/N2d)
            Rinv_actual = (ftool.power_1d3[1:].sum()/Nz)/(ftool.power_1d2[1:].sum()/N2d)
            sigma_Brunt = (sigma_col*Rinv).real
            sigma_Brunt_actual = (sigma_col*Rinv_actual).real
            #print("Brunt", sigma_Brunt/sigma_rho)
            #print("One", sigma_Brunt_actual/sigma_rho)

            R1 =sigma_Brunt/sigma_rho
            R2 = sigma_Brunt_actual/sigma_rho
            #ax.set(title='%s %0.1e %0.1e'%(sim,R2,R1))
            ax.text(0.1,0.1,"%s"%sim,transform=ax.transAxes)
            ax.text(0.25,0.1,"%0.1f"%(np.sqrt(R1)),transform=ax.transAxes)

        if 1:
            xvals = this_sim.avg_spectra['k2d']
            mask = slice(4,25)

            N2d = ftool.rho2.size
            mean_column = ftool.rho2.sum()/N2d
            sigma_col = ((ftool.rho2-mean_column)**2).sum().real/N2d
            sigma_2d_fft = ftool.power_1d2[mask].sum().real/N2d
            #print("col %0.2e fft %0.2e ratio %0.2e"%(sigma_col, sigma_2d_fft, sigma_col/sigma_2d_fft))
            Rinv = (( ftool.k2d*ftool.power_1d2)[mask].sum()/Nz)/((ftool.power_1d2)[mask].sum()/N2d)
            Rinv_actual = (ftool.power_1d3[mask].sum()/Nz)/(ftool.power_1d2[mask].sum()/N2d)
            sigma_Brunt = (sigma_col*Rinv).real
            sigma_Brunt_actual = (sigma_col*Rinv_actual).real
            #print("Brunt", sigma_Brunt/sigma_rho)
            #print("One", sigma_Brunt_actual/sigma_rho)

            R1 =sigma_Brunt/sigma_rho
            R2 = sigma_Brunt_actual/sigma_rho
            #ax.set(title='%s %0.1e %0.1e'%(sim,R2,R1))
            ax.text(0.35,0.1,"%0.1f"%(np.sqrt(R1)),transform=ax.transAxes)
            ax.text(0.45,0.1,"%0.1f"%(np.sqrt(sigma_rho)),transform=ax.transAxes)
            ax.text(0.55,0.1,"%0.3f"%(np.sqrt(Rinv.real)),transform=ax.transAxes)

        if 1:
            a2.scatter(np.sqrt(Rinv.real), np.sqrt(sigma_rho), c=[this_sim.color],s=this_sim.marker_size*100)
            a2b.scatter(this_sim.Ms_mean, np.sqrt(sigma_rho), c=[this_sim.color],s=this_sim.marker_size*100)
            SR,SB=np.sqrt(sigma_rho), np.sqrt(sigma_Brunt)
            a2c.scatter(SR,1-SB/SR ,c=[this_sim.color],s=this_sim.marker_size*100)



        
        #a2c.plot([0,2.5],[0,2.5],c=[0.5]*3)
        f2.savefig('%s/Rsig'%dl.plotdir)
        print(outname)
        #ax.legend(loc=0)
        #ax.set(yscale='log',xscale='log')

    outname="%s/brunt_spectra"%dl.plotdir
    print(outname)

    fig.tight_layout()
    fig.savefig(outname)
    return ftool

def drill(simlist, ftool=None):
    plt.close('all')

    nrow = int(np.ceil(len(simlist)/3))
    ncol = min([3,len(simlist)])
    fig,axes = plt.subplots(nrow,ncol,figsize=(12,12))
    axlist=axes.flatten()
    #fig.subplots_adjust(wspace=0, hspace=0)

    for i,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]
        this_sim.load()
        simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)                                                                                                                 #ax.scatter( raq.quan3[sim]['maavg'], raq.quan3[sim]['msavg'],c=sim_colors.color[sim], marker=sim_colors.marker[sim],s=60)
        frame=this_sim.ann_frames[-1]


        
        ftool=bucket.things.get(sim,None)
        if ftool is None:
            rho_full, rho = bt.get_cubes(sim,frame)

            ftool = bt.fft_tool(rho)
            ftool.do3()
            ftool.do2(projax=0)
        bucket.things[sim]=ftool

        ax=axlist[i]
        ax.plot(ftool.k2d, ftool.power_1d2.real,c='r', label='P2d')
        ax.plot(ftool.k3d, ftool.power_1d3.real,c='g', label='P3d')
        ax.plot(ftool.k2d,ftool.k2d*ftool.power_1d2.real,c='b', label = 'k P2d')
