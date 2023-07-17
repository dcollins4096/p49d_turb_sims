
from GL import *

import simulation
def plot_all_spectra(simlist):
    for this_simname in simlist:
        print('plot',this_simname)
        this_sim = simulation.corral[this_simname]
        this_sim.read_all_spectra()
        this_sim.fit_all_spectra()

        #
        # Primitive 
        #

        frames = this_sim.all_frames
        
        fig,axes=plt.subplots(7,3, figsize=(8,12))
        #ax0=axes[0];ax1=axes[1];ax2=axes[2]

        avg_slope={}
        for field in ['density','velocity','Htotal']:
            avg_slope[field]=0
        for axis in 'xyz':
            avg_slope[axis]={}
            for field in ['ClTT','ClEE','ClBB']:
                avg_slope[axis][field]=0
        for field in ['density','velocity','Htotal']:

            NF = 0
            for frame in frames:
                this_slope=this_sim.slopes[frame][field]
                if this_slope is None:
                    continue
                avg_slope[field]+=this_slope
                NF+=1
            avg_slope[field]/=NF
        for axis in 'xyz':
            for field in ['ClTT','ClEE','ClBB']:
                for frame in frames:
                    this_slope=this_sim.slopes[frame][axis][field]
                    if this_slope is None:
                        continue
                    avg_slope[axis][field]+=this_slope
                    NF+=1
                avg_slope[axis][field]/=NF

        for frame in frames:
            aaa = this_sim.all_spectra[frame]
            k = aaa['k3d']
            for nf,field in enumerate(['density','velocity','Htotal']):
                thax=axes[0][nf]
                slope = avg_slope[field]
                spec=aaa[field]
                ok = (aaa[field]>0)*(k>0)
                comp=k[ok]**-slope
                #comp=1
                thax.plot( k[ok], spec[ok]*comp, c=[0.5]*4)
            for nax,axis in enumerate('xyz'):
                k2d = aaa[axis]['k2d']
                for nf,field in enumerate(['ClTT','ClEE','ClBB','ClTE','ClTB','ClEE']):
                    thax=axes[nf+1][nax]
                    ok = k2d>0
                    spec=aaa[axis][field]
                    thax.plot( k2d[ok], spec[ok])


        for nf,field in enumerate(['density','velocity','Htotal']):
            thax=axes[0][nf]
            thax.set(xscale='log',yscale='log',xlabel='K',ylabel='%s Power'%field, 
                title="%0.2f"%avg_slope[field])
            for nax,axis in enumerate('xyz'):
                for nf,field in enumerate(['ClTT','ClEE','ClBB']):
                    thax=axes[nf+1][nax]
                    thax.set(xscale='log',yscale='log',xlabel='k2d',ylabel='%s %s'%(field,axis),
                             title="%0.2f"%avg_slope[axis][field])
                for nf,field in enumerate(['ClTE','ClTB','ClEB']):
                    thax=axes[nf+4][nax]
                    thax.set(xscale='log',xlabel='k2d',ylabel='%s %s'%(field,axis))
                    thax.set_yscale('symlog',linthresh=1e-2)

        fig.tight_layout()
        fig.savefig('%s/all_spectra_%s'%(dl.plotdir,this_sim.name))




