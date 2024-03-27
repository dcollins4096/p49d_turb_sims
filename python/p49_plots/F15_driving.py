
from GL import *

plotdir = "%s/PigPen"%os.environ['HOME']
import simulation

def plot_accel(simlist):
    for ns,sim in enumerate(simlist):
        plt.close('all')
        fig,axes=plt.subplots(1,1)
        ax0=axes
        this_sim = simulation.corral[sim]
        this_sim.load()
        

        xvals = this_sim.avg_spectra['k3d']

        for frame in this_sim.ann_frames:
            yvals=this_sim.all_spectra[frame]['acceleration']
            try:
                len(yvals)
                ax0.plot(xvals, yvals,c=this_sim.color,linestyle=this_sim.linestyle)
            except:
                continue
        ax0.set(xscale='log',yscale='log')
        fig.savefig('%s/driving_%s.png'%(plotdir,sim))

def plot_driving(simlist):

    total_columns = min([2,len(simlist)])
    total_rows = len(simlist)//total_columns
    fig1,axes1=plt.subplots(total_rows, total_columns, figsize=(20,20))
    fig2,axes2=plt.subplots(total_rows, total_columns, figsize=(20,20))
    for nsim,sim in enumerate(simlist):
        nr = nsim//total_rows
        nc = nsim % total_columns
        this_sim = simulation.corral[sim]
        ax1 = axes1[nr][nc]
        ax2 = axes2[nr][nc]
        #ax1=axes1
        #ax2=axes2
        aaa=[ax1,ax2]

        frame0 = this_sim.ann_frames[0]
        frame1 = this_sim.ann_frames[-1]
        proj_ax=0
        field = 'x-acceleration'
        for nf, frame in enumerate([frame0,frame1]):
            ds = this_sim.load_ds(frame)
            proj = ds.proj(field,proj_ax)
            #pw = proj.to_pw()
            frb = proj.to_frb(1,[512,512])
            aaa[nf].imshow(frb[field])
    fig1.savefig('%s/drive1'%plotdir)
    fig2.savefig('%s/drive2'%plotdir)




