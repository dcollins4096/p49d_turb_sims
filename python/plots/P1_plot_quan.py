from GL import *
reload(sim_colors)
import data_locations as dl


plotdir = dl.plotdir

import read_avg_quan as raq
#reload(raq)

"""
alf_x_avg                Dataset {1}
alf_x_std                Dataset {1}
alf_y_avg                Dataset {1}
alf_y_std                Dataset {1}
alf_z_avg                Dataset {1}
alf_z_std                Dataset {1}
bx_avg                   Dataset {1}
bx_std                   Dataset {1}
by_avg                   Dataset {1}
by_std                   Dataset {1}
bz_avg                   Dataset {1}
bz_std                   Dataset {1}
density_avg              Dataset {1}
density_std              Dataset {1}
time                     Dataset {1}
vx_avg                   Dataset {1}
vx_std                   Dataset {1}
vy_avg                   Dataset {1}
vy_std                   Dataset {1}
vz_avg                   Dataset {1}
vz_std                   Dataset {1}
"""
reload(raq)
#sim_list=sim_colors.simlist
#sim_list=['6_1']
import simulation as sim
def plot_quan(sim_list):
    plt.close('all')
    fig,ax=plt.subplots(3,4,figsize=(12,8))
    for ns,sim_name in enumerate(sim_list):
        raq.read(sim_name)
        this_sim=sim.corral[sim_name]
        time = raq.quan_time[sim_name]['time']+0
        print("%10s max %0.2f tdyn %0.3f t/tdyn %0.3f"%(sim_name,time.max(), this_sim.tdyn, time.max()/ this_sim.tdyn))
        time /= this_sim.tdyn
        print(time.max())
        #time = nar(range(len(raq.quan_time[sim]['time'])))
        #print(time)
        QQQ = raq.quan_time[sim_name]
        vx_avg = QQQ['vx_avg']
        vy_avg = QQQ['vy_avg']
        vz_avg = QQQ['vz_avg']
        vmag = (vx_avg**2+vy_avg**2+vz_avg**2)
        ax[0][0].plot( time, QQQ['vx_avg'], c=this_sim.color)
        ax[0][1].plot( time, QQQ['vy_avg'], c=this_sim.color)
        ax[0][2].plot( time, QQQ['vz_avg'], c=this_sim.color)
        ax[0][3].plot( time, vmag, c=this_sim.color)
        if 0:
            fig2,ax2=plt.subplots(3,4,figsize=(12,8))
            ax2[0][0].plot( time, QQQ['vx_avg'], c=sim_colors.color[sim])
            ax2[0][1].plot( time, QQQ['vy_avg'], c=sim_colors.color[sim])
            ax2[0][2].plot( time, QQQ['vz_avg'], c=sim_colors.color[sim])
            ax2[0][3].plot( time, vmag, c=sim_colors.color[sim])
            fig2.savefig('%s/quan_mon_%s.png'%(plotdir,sim))
            plt.close(fig2)


    fig.savefig('%s/quan_monster.pdf'%plotdir)

