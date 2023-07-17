from GL import *
reload(sim_colors)
import data_locations as dl


plotdir = dl.plotdir

#import read_avg_quan as raq
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
#reload(raq)
#sim_list=sim_colors.simlist
#sim_list=['6_1']
import simulation as sim
def plot_quan(sim_list):
    plt.close('all')
    fig,ax=plt.subplots(3,4,figsize=(12,8))
    if len(sim_list)>1:
        outname = '%s/avg_quan_multi'%(dl.plotdir)
    else:
        outname = '%s/avg_quan_%s'%(dl.plotdir,sim_list[0])
    for ns,sim_name in enumerate(sim_list):
        this_sim=sim.corral[sim_name]
        this_sim.read_avg_quan()

        time = this_sim.quan_time['time']+0
        print("%10s max %0.2f tdyn %0.3f t/tdyn %0.3f"%(sim_name,time.max(), this_sim.tdyn, time.max()/ this_sim.tdyn))
        time /= this_sim.tdyn
        #time = nar(range(len(raq.quan_time[sim]['time'])))
        #print(time)
        QQQ = this_sim.quan_time
        vx_avg = QQQ['vx_avg']
        vy_avg = QQQ['vy_avg']
        vz_avg = QQQ['vz_avg']
        vmag = (vx_avg**2+vy_avg**2+vz_avg**2)
        ax[0][0].plot( time, QQQ['vx_avg'], c=this_sim.color)
        ax[0][1].plot( time, QQQ['vy_avg'], c=this_sim.color)
        ax[0][2].plot( time, QQQ['vz_avg'], c=this_sim.color)
        ax[0][3].plot( time, vmag, c=this_sim.color)
        ax[0][0].set(xlabel='t/tdyn',ylabel=r'$\langle v_x \rangle$')
        ax[0][1].set(xlabel='t/tdyn',ylabel=r'$\langle v_y \rangle$')
        ax[0][2].set(xlabel='t/tdyn',ylabel=r'$\langle v_z \rangle$')
        ax[0][3].set(xlabel='t/tdyn',ylabel=r'||$\langle v_i \rangle$||')

        bx_avg = QQQ['bx_avg']
        by_avg = QQQ['by_avg']
        bz_avg = QQQ['bz_avg']
        bmag = (bx_avg**2+by_avg**2+bz_avg**2)
        ax[1][0].plot( time, QQQ['bx_avg'], c=this_sim.color)
        ax[1][1].plot( time, QQQ['by_avg'], c=this_sim.color)
        ax[1][2].plot( time, QQQ['bz_avg'], c=this_sim.color)
        ax[1][3].plot( time, bmag, c=this_sim.color)
        ax[1][0].set(xlabel='t/tdyn',ylabel=r'$\langle b_x \rangle$')
        ax[1][1].set(xlabel='t/tdyn',ylabel=r'$\langle b_y \rangle$')
        ax[1][2].set(xlabel='t/tdyn',ylabel=r'$\langle b_z \rangle$')
        ax[1][3].set(xlabel='t/tdyn',ylabel=r'||$\langle b_i \rangle$||')

        ax[2][0].plot(time,QQQ['vrms'], c=this_sim.color)
        ax[2][1].plot(time,QQQ['ma'], c=this_sim.color)
        ms = this_sim.quan3['msavg']
        ax[2][0].axhline(ms, label="%0.1f"%ms)
        ax[2][0].set(ylabel='vrms')
        ax[2][0].legend(loc=0)
        ax[2][1].set(ylabel=r'$v_{rms}/\langle B \rangle/"\sqrt{4\pi}"')
        ax[2][1].axhline(this_sim.quan3['maavg'], label="%0.1f"%this_sim.quan3['maavg'])
        ax[2][1].legend(loc=0)

    fig.tight_layout()
    fig.savefig(outname)
    print(outname)

