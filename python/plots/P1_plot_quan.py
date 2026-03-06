from GL import *
reload(sim_colors)
import data_locations as dl


plotdir = dl.plot_dir

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
def plot_all_mach(sim_list, ncol=4):
    plt.close('all')
    nrow = max( len(sim_list)//ncol, 1)
    fig,ax=plt.subplots(nrow,ncol,figsize=(12,8))
    fig.subplots_adjust(wspace=0, hspace=0)
    if len(sim_list)>1:
        outname = '%s/avg_quan_multi.pdf'%(plotdir)
    else:
        outname = '%s/avg_quan_%s.pdf'%(plotdir,sim_list[0])
    for ns,sim_name in enumerate(sim_list):
        nx = ns//ncol
        ny = ns%ncol
        this_sim=sim.corral[sim_name]
        this_sim.read_avg_quan()

        time = this_sim.quan_time['time']+0
        print("%10s max %0.2f tdyn %0.3f t/tdyn %0.3f"%(sim_name,time.max(), this_sim.tdyn, time.max()/ this_sim.tdyn))
        time /= this_sim.tdyn
        #time = nar(range(len(raq.quan_time[sim]['time'])))
        #print(time)
        QQQ = this_sim.quan_time
        #vx_avg = QQQ['vx_avg']
        #vy_avg = QQQ['vy_avg']
        #vz_avg = QQQ['vz_avg']
        #vmag = (vx_avg**2+vy_avg**2+vz_avg**2)
        ax[nx][ny].plot(time,QQQ['vrms']/np.sqrt(3), c=this_sim.color)
        ms = this_sim.quan_mean['msavg']/np.sqrt(3)
        ax[nx][ny].axhline(ms, label="avg = %0.1f"%ms)
        ax[nx][ny].axhline(this_sim.Ms_nom, label='Nominal = %0.1f'%this_sim.Ms_nom, c='r')
        #ax[nx][ny].legend(loc=1)
        #ax[nx][ny].set(xlabel='t/tdyn', ylabel='Mach')
        if ny == 0:
            ax[nx][ny].set(ylabel='1d Mach')
        else:
            ax[nx][ny].set(yticks=[])
        if nx == nrow-1:
            ax[nx][ny].set(xlabel='t/tdyn')
        else:
            ax[nx][ny].set(xticks=[])
        #ax[nx][ny].set(ylim=[0,10])

    fig.savefig('%s/all_mach'%plot_dir)

def plot_quan(sim_list):
    plt.close('all')
    fig,ax=plt.subplots(2,4,figsize=(12,8))
    if len(sim_list)>1:
        outname = '%s/avg_quan_multi.pdf'%(plotdir)
    else:
        outname = '%s/avg_quan_%s.pdf'%(plotdir,sim_list[0])
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
        #ax[0][3].plot( time, vmag, c=this_sim.color)
        ax[0][0].set(xlabel='t/tdyn',ylabel=r'$\langle v_x \rangle$')
        ax[0][1].set(xlabel='t/tdyn',ylabel=r'$\langle v_y \rangle$')
        ax[0][2].set(xlabel='t/tdyn',ylabel=r'$\langle v_z \rangle$')
        ax[0][3].set(xlabel='t/tdyn',ylabel=r'||$\langle v_i \rangle$||')

        ax[0][3].plot(time,QQQ['vrms'], c=this_sim.color)
        ms = this_sim.quan_mean['msavg']
        ax[0][3].axhline(ms, label="avg = %0.1f"%ms)
        ax[0][3].axhline(this_sim.Ms_nom, label='Nominal = %0.1f'%this_sim.Ms_nom, c='r')
        ax[0][3].set(ylabel='vrms')
        ax[0][3].legend(loc=0)

        if this_sim.do_magnetic:
            bx_avg = QQQ['bx_avg']
            by_avg = QQQ['by_avg']
            bz_avg = QQQ['bz_avg']
            #bmag = (bx_avg**2+by_avg**2+bz_avg**2)
            ax[1][0].plot( time, QQQ['bx_avg'], c=this_sim.color)
            ax[1][0].axhline( this_sim.B_nom, label='Nominal = %0.2f'%this_sim.B_nom, c='r')
            ax[1][1].plot( time, QQQ['by_avg'], c=this_sim.color)
            ax[1][2].plot( time, QQQ['bz_avg'], c=this_sim.color)
            #ax[1][3].plot( time, bmag, c=this_sim.color)
            ax[1][0].set(xlabel='t/tdyn',ylabel=r'$\langle b_x \rangle$')
            ax[1][1].set(xlabel='t/tdyn',ylabel=r'$\langle b_y \rangle$')
            ax[1][2].set(xlabel='t/tdyn',ylabel=r'$\langle b_z \rangle$')

            ax[1][3].plot(time,QQQ['ma'], c=this_sim.color)
            ax[1][3].set(ylabel=r'$v_{rms}/\langle B \rangle/\sqrt{4\pi}$')
            ax[1][3].axhline(this_sim.quan_mean['maavg'], label="avg=%0.1f"%this_sim.quan_mean['maavg'])
            ax[1][3].axhline(this_sim.Ma_nom, label = "Nominal = %0.1f"%this_sim.Ma_nom, c='r')
            ax[1][3].legend(loc=0)

    fig.tight_layout()
    fig.savefig(outname)
    print(outname)

