
from GL import *
from downsample import volavg



import simulation
reload(simulation)
import simulation_info.all_sims
#import old_brunt_tools as bt
import dtools.math.brunt_tools as bt
reload(bt)

import p68_laser.saver as saver
def plot_all_brunt(sim_list, projax=0):

    ncol = 3
    nrow = np.ceil(len(sim_list)/3).astype('int')
    fig,axes = plt.subplots(nrow,ncol,figsize=(12,12))

    if nrow==1:
        axes=[axes]
    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        frame = this_sim.ann_frames[-1]
        if sim not in saver.bucket:
            rho = this_sim.load_small_rho(frame)
            ftool = bt.fft_tool(rho)
            ftool.do3()
            ftool.do2(projax=projax)
            saver.bucket[sim]=ftool
        else:
            ftool=saver.bucket[sim]
        nc = nsim%ncol
        nr = nsim//ncol
        ax = axes[nr][nc]
        bt.plot_brunt(ftool,method='full',ax=ax)
    for ax in axes.flatten():
        ax.set(xticks=[],yticks=[])
    fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,top=1,bottom=0)
    #fig.tight_layout()
    fig.savefig('%s/all_brunt'%(plotdir))


def plot_sigmas(sim_list, projax=0):

    fig,axes = plt.subplots(2,2)

    for nsim,sim in enumerate(sim_list):
        this_sim=simulation.corral[sim]
        this_sim.load()
        frame = this_sim.ann_frames[-1]
        if sim not in saver.bucket:
            rho = this_sim.load_small_rho(frame)
            ftool = bt.fft_tool(rho)
            ftool.do3()
            ftool.do2(projax=projax)
            saver.bucket[sim]=ftool
        else:
            ftool=saver.bucket[sim]
        axes[0][0].scatter(this_sim.Ms_mean,1-ftool.sigma_x3d/ftool.sigma_k3d, c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[0][0].set(xlabel='Ms',ylabel=r'$1-\sigma_{3x}/\sigma_{3k}$', title='3x vs 3k')
        axes[0][1].scatter(this_sim.Ms_mean,1-ftool.sigma_x2d/ftool.sigma_k2d, c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[0][1].set(xlabel='Ms',ylabel=r'$1-\sigma_{2x}/\sigma_{2k}$', title='2x vs 2k')
        axes[1][0].scatter(this_sim.Ma_mean, 1-ftool.sigma_k3d/ftool.sigma_k2dk, c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[1][0].set(xlabel=r'$M_a$',ylabel=r'$1-\sigma_{3k}/\sigma_{k2k}$', title='3k vs k2k')
        axes[1][1].scatter( this_sim.Ma_mean, ftool.ratio_1,c=[this_sim.color], marker=this_sim.marker,s=this_sim.marker_size*100)
        axes[1][1].set(xlabel=r'$M_a$',ylabel='ratio',title='goal')


    fig.tight_layout()
    fig.savefig('%s/all_sigmas'%(plotdir))


