
from GL import *

import simulation
def plot_all_spectra(simlist, all_or_ann='all'):
    """ all frames or ann (analysis) frames"""
    plt.close('all')
    for this_simname in simlist:
        print('plot',this_simname)
        this_sim = simulation.corral[this_simname]
        this_sim.read_all_spectra()
        this_sim.fit_all_spectra()

        #
        # Primitive 
        #

        if all_or_ann == 'all':
            frames = this_sim.all_frames
        elif all_or_ann == 'ann':
            frames = this_sim.ann_frames
        else:
            print("ERROR. all_or_ann must be all or ann.")
            pdb.set_trace()
        
        fig,axes=plt.subplots(7,3, figsize=(8,12))
        #ax0=axes[0];ax1=axes[1];ax2=axes[2]

        avg_slope={}
        for field in this_sim.products_positive:
            avg_slope[field]=0
            NF = 0
            for frame in frames:
                this_slope=this_sim.slopes[frame][field]
                if this_slope is None:
                    continue
                avg_slope[field]+=this_slope
                NF+=1
            avg_slope[field]/=NF

        for frame in frames:
            aaa = this_sim.all_spectra[frame]
            for nf,field in enumerate(this_sim.products):
                if field in this_sim.products_3d:
                    xvals = aaa['k3d']
                else:
                    xvals = aaa['k2d']
                ncol = nf%3
                nrow = nf//3
                thax=axes[nrow][ncol]
                spec=aaa[field]
                #ok = (aaa[field]>0)*(xvals>0)
                ok = (xvals>0)
                if field in this_sim.products_positive:
                    slope = avg_slope[field]
                    spec[ok] = spec[ok]*xvals[ok]**np.abs(slope)
                #comp=1
                thax.plot( xvals[ok], spec[ok], c=[0.5]*4)


        for nf,field in enumerate(this_sim.products):
            ncol = nf%3
            nrow = nf//3
            thax=axes[nrow][ncol]
            if field in this_sim.products_3d:
                thax.set(xscale='log',yscale='log',xlabel='K',ylabel='%s Power'%field, 
                    title="%0.2f"%avg_slope[field])
            elif field in this_sim.products_positive:
                thax.set(xscale='log',yscale='log',xlabel='k2d',ylabel='%s'%(field),
                         title="%0.2f"%avg_slope[field])
            else:
                thax.set(xscale='log',xlabel='k2d',ylabel='%s'%(field))
                thax.set_yscale('symlog',linthresh=1e-2)

        fig.tight_layout()
        fig.savefig('%s/all_spectra_%s'%(dl.plotdir,this_sim.name))




