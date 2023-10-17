from GL import *
import scipy


import simulation


def plot_dt(sim_list):

    fig,axes=plt.subplots(1,3)
    ax=axes[0];ax1=axes[1];ax2=axes[2]
    va = []
    vs = []
    dts = []
    tmax=[]
    nmax=[]
    dwall_all=[]
    meanwall=[]
    totalwall=[]
    tmin=[]
    for ns,sim in enumerate(sim_list):
        this_sim = simulation.corral[sim]
        this_sim.load()
        output_log = "%s/OutputLog"%this_sim.data_location
        fptr = open(output_log)
        lines=fptr.readlines()
        fptr.close()

        cycle=[]
        time=[]
        frames=[]
        wall=[]

        for line in lines:
            lll = line.split()
            frame= int( lll[2][-4:])
            if frame not in this_sim.ann_frames:
                continue
            if frame in frames:
                continue
            frames.append(frame)
            cycle.append(int(lll[3]))
            time.append(float(lll[4]))
            wall.append(float(lll[5]))

        cycle=nar(cycle)
        time=nar(time)
        wall=nar(wall)

        dt = (time[1:]-time[:-1])
        dc = (cycle[1:]-cycle[:-1])
        dwall=(wall[1:]-wall[:-1])
        dwall_dc = dwall/dc
        prob_dwall = dwall_dc+0
        prob_dwall.sort()
        mean_dwall = prob_dwall[:prob_dwall.size//2].mean()
        #mean_dwall = prob_dwall.mean()
        meanwall.append(mean_dwall)
        total_wall = dwall[ dwall_dc < 20*mean_dwall].sum()
        totalwall.append(total_wall)

        print(this_sim.name, mean_dwall, total_wall/3600)
        #print( (dwall_dc>=20*mean_dwall).sum())

        dwall_all+=list(dwall)
        dtdc=dt/dc
        tc = 0.5*(time[1:]+time[:-1])
        ax.plot(dtdc,c=this_sim.color)

        ms = this_sim.quan_mean['msavg']
        ma = this_sim.quan_mean['maavg']
        dtmean = dtdc.mean()
        vs.append(ms)
        va.append(ma)
        dts.append(dtmean)
        tmax.append(max(time))
        tmin.append(min(time))
        nmax.append(max(cycle))
        #ax1.scatter(ms,dtmean,color=this_sim.color)

        #lastframe=max(frames)
        #ds = yt.load("%s/DD%04d/data%04d"%(this_sim.data_location,lastframe,lastframe))
        #print(ds['CourantSafetyNumber'])
        #ax2.hist( np.log10( dwall/dwall.mean()), histtype='step')
        if 0:
            x = dwall+0
            x.sort()
            y=(np.arange(x.size)/x.size)[::-1]
            #ax2.plot(x,y)
            ax2.plot(x[::-1])
            ax2.set(xlim=[-1,10], yscale='log')
        if 0:
            dn = np.arange(dwall.size)
            ind = np.argmax(dwall)
            ax2.plot(dn-ind, dwall)
            gotit=True
            ax2.set(xscale='log')

    #ax2.hist(meanwall)
    meanwall=nar(meanwall)
    ncore=4096
    zucs=(512**3)/meanwall/ncore
    ax2.hist(zucs)
    vs=nar(vs)
    va=nar(va)
    dts=nar(dts)
    dx=1./512
    dt2 = 0.5*dx/(1+vs)
    #ax1.plot( vs, dt2, c='k')
    ytmp=np.log10(dt2)
    ax1.scatter( vs, dts)
    ax1.set(yscale='log')
    def dtthing(x,m,b):
        #return np.exp(m*x+b)
        return (m*x+b)
    def dtthinge(x,m,b):
        return np.exp(m*x+b)
        #return (m*x+b)
    B=vs/va
    if 1:
        pfit, pcov = scipy.optimize.curve_fit( dtthing, vs, np.log(dts))
        print(pfit)
        ax1.plot( vs, np.exp(dtthing( vs, *pfit)), c='r')
        ax1.plot( vs, 0.1*dx/(1+vs+B),c='m')
        print(pfit[0])
        #ax1.plot( vs,  pfit[0]*vs+pfit[1], c='r')
        fig.savefig("plots_to_sort/dt")

    args = np.argsort(dts)
    for i in args:
        est_wall = totalwall[i]*(tmax[i]/(tmax[i]-tmin[i]))
        #print("mach %0.2f alf %0.2f B %0.2f dt %0.2e Tmax %0.2f Nsteps %0.2e Twall %0.2f est %0.2f"%(vs[i], va[i], B[i], dts[i], tmax[i], nmax[i], totalwall[i]/3600, est_wall/3600))
        zucs=512**3*nmax[i]/totalwall[i]/4096
        print("mach %0.2f alf %0.2f Tmax %0.2f Nsteps %0.2e Twall %0.2f est %0.2f"%(vs[i], va[i], tmax[i], nmax[i], totalwall[i]/3600, zucs))


