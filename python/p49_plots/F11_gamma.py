from GL import *

import simulation

def plot_gamma(simlist):

    plt.close('all')
    fig,axes=plt.subplots(2,2)
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]
    all_gamma=[]

    fig2,aaa2=plt.subplots(1,1)
    fig3,aaa3=plt.subplots(1,1)
    fig4,aaa4=plt.subplots(1,2)
    
    for ns,sim in enumerate(simlist):
        this_sim = simulation.corral[sim]
        this_sim.load()
        simulation.set_colors(this_sim,cmap_name='gist_rainbow')
        xvals = this_sim.avg_spectra['k2d']
        EE= this_sim.avg_spectra['ClEEy'] 
        BB= this_sim.avg_spectra['ClBBy']
        EB= this_sim.avg_spectra['ClEBy']
        fit_range =this_sim.get_fitrange(xvals)
        mask = (xvals > fit_range[0])*(xvals < fit_range[1])

        kwargs={'color':this_sim.color,'linestyle':this_sim.linestyle}
        ax0.plot( xvals[mask],EE[mask],**kwargs)
        ax1.plot( xvals[mask], (EE-BB)[mask], **kwargs)
        ax2.plot( xvals[mask], EB[mask],**kwargs)

        Top,Bot= 2*EB[mask].mean(),(EE-BB)[mask].mean()
        gamma_m = Top/Bot
        #my_gamma = np.arcsin(gamma_m)/4#*180/(np.pi)
        #my_gamma = np.arcsin(gamma_m)/4
        my_gamma = gamma_m
        #aaa4[0].scatter(this_sim.Ms_mean, my_gamma, **kwargs)
        aaa4[0].scatter(this_sim.Ms_mean, my_gamma, **kwargs)
        aaa4[0].set(xlabel='Ms',ylabel='mean gamma')
        if 1:
            #aaa4[1].scatter( my_gamma,Top,**kwargs)
            aaa4[1].scatter( Top,1/Bot,**kwargs)
            aaa4[1].set(xlabel='EB',ylabel='1/EE-BB',yscale='log', xscale='log')
        #aaa4[1].scatter( gamma_m, my_gamma, **kwargs)

        if my_gamma>0.2:
            print(my_gamma)
        eb_fid = np.geomspace(1e-3,0.5,64)
        y = gamma_m/eb_fid
        aaa4[1].plot(eb_fid,y, **kwargs)



        gamma1= 2*EB/(EE-BB)
        newmask = (np.abs(gamma1)<=1)*mask
        gamma = np.arcsin(gamma1[newmask])/4
        gamma *= 180/np.pi
        all_gamma+=list(gamma)
        #pdb.set_trace()
        #ax3.plot( xvals[newmask],gamma,**kwargs)
        #kwargs['histtype']='step'
        #kwargs['bins']=np.linspace(-0.2,0.2,16)
        #kwargs['bins']=np.geomspace(1e-3,0.2,16)

        #ax3.hist(np.abs(gamma),**kwargs)
        #ax3.hist(gamma,**kwargs)
        #aaa2.hist(gamma, **kwargs)

        gp = gamma[gamma>=0]
        gp.sort()
        y = (np.arange(gp.size)/gp.size)[::-1]
        aaa2.plot(gp,y,**kwargs)
        gp = gamma[gamma<0]
        gp.sort()
        y = (np.arange(gp.size)/gp.size)

        aaa2.plot(gp,y,**kwargs)

        kwargs['s']=1
        #aaa3.scatter(gamma,180/(np.pi*4)* gamma1[newmask], **kwargs)
        f=180/(4*np.pi)
        #aaa3.scatter(gamma,f* gamma1[newmask], **kwargs)
        #aaa3.scatter(gamma, f*(2*EB/(EE-BB))[newmask],**kwargs)
        #aaa3.scatter(gamma, ((EE-BB))[newmask],**kwargs)
        #aaa3.set(yscale='log')
        #aaa3.scatter(EB[newmask], gamma,**kwargs)
        aaa3.scatter(EB[newmask], (1/(EE-BB))[newmask],**kwargs)
        #aaa3.scatter(EB[newmask], ((EE-BB))[newmask],**kwargs)

    all_gamma=nar(all_gamma)
    #all_gamma*=180/np.pi
    #minmin=2e-3
    #bins = np.geomspace(minmin,0.2,16)
    #bins = np.concatenate([-bins[::-1],bins])
    #ax3.hist(all_gamma, histtype='step',bins=bins)
    stuffs = [all_gamma[all_gamma>=0], np.abs(all_gamma[all_gamma<0])]
    for n,G in enumerate(stuffs):
        #G = all_gamma+0
        G.sort()
        y = np.arange(G.size)/G.size
        #ax3.plot(G,y, c='rg'[n])

    #aaa2.hist(all_gamma, histtype='step',color='k')
    fig2.savefig('plots_to_sort/gamma_hist')
    fig3.savefig('plots_to_sort/gamma_play')
    fig4.savefig('plots_to_sort/gamma_mean')



    ax0.set(xscale='log',yscale='log', title='EE')
    ax1.set(xscale='log',yscale='log', title='EE-BB')
    ax2.set(xscale='log', title='EB')
    ax2.set_yscale('symlog',linthresh=0.1)
    #ax3.set(xscale='log')
    #ax3.set_yscale('symlog',linthresh=0.1)
    ax3.set(xscale='log')
    #ax3.set_xscale('symlog',linthresh=1)


    fig.savefig('plots_to_sort/gamma')
