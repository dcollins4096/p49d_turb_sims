
from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(sim_colors)
reload(queb3)
import get_all_quantities as gaq
reload(gaq)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)


#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
#reload(read_stuff)  
spectra_dict = rs.spectra_dict
quan3 = rs.quan3



#
# Spectra plots
#

plotdir=os.environ['HOME']+"/PigPen"
from scipy.stats import ks_2samp
def do_ks_string(test_dist, original_dist):
        KS_output = ks_2samp(test_dist, original_dist)
        crit_stat = 1.36*np.sqrt( (original_dist.size + test_dist.size)/(original_dist.size*test_dist.size))
        ks_string = r'$D=%0.2f (D_c=%0.2f, p=%f)$'%(KS_output.statistic,crit_stat, KS_output.pvalue)
        return ks_string
def do_ks(test_dist, original_dist):
        KS_output = ks_2samp(test_dist, original_dist)
        crit_stat = 1.36*np.sqrt( (original_dist.size + test_dist.size)/(original_dist.size*test_dist.size))
        #ks_string = r'$D=%0.2f (D_c=%0.2f, p=%f)$'%(KS_output.statistic,crit_stat, KS_output.pvalue)
        return KS_output.statistic,crit_stat, KS_output.pvalue

class slope_bucket():
    def __init__(self):
        self.spext={}
        self.speyt={}
        self.slope_array={}
        self.chi_array={}
        self.slopext=defaultdict(dt.extents)
        self.chiext=defaultdict(dt.extents)
        self.amp16={}
        self.k16={}

field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
if 0:
    #
    # Development code.  Now part of read_stuff
    #
    #fill objects with slope array
    slope_master = dt.extents()
    spec_master = dt.extents()
    for ns,sim in enumerate(sim_colors.simlist):
        spectra_dict[LOS][sim].slope_bucket = slope_bucket()
    for nf,field in enumerate(field_list):
        for ns,sim in enumerate(sim_colors.simlist):
            print(field,sim)
            for LOS in 'xyz':
                proj=rs.proj_dict[LOS][sim]
                x0_list = proj.lcent[2:10:2]
                x1_list = proj.lcent[11:51:4]
                max_npoints=np.argmax(x1_list)-np.argmin(x0_list)+1
                spectra_dict[LOS][sim].slope_bucket.spext[field]=dt.extents()
                spectra_dict[LOS][sim].slope_bucket.speyt[field]=dt.extents()
                spectra_dict[LOS][sim].slope_bucket.slope_array[field]= np.zeros([x0_list.size, x1_list.size]) +10
                spectra_dict[LOS][sim].slope_bucket.chi_array[field]= np.zeros([x0_list.size, x1_list.size]) +10

                spec=spectra_dict[LOS][sim].spectra[field]
                spectra_dict[LOS][sim].slope_bucket.amp16[field] = spec[16]
                spectra_dict[LOS][sim].slope_bucket.k16[field] = proj.lcent[16]
                for n0,x0 in enumerate(x0_list):
                    for n1,x1 in enumerate(x1_list):

                        ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                        xvals = proj.lcent[ok]
                            
                        spec=spectra_dict[LOS][sim].spectra[field][ok]
                        spec_master(spec)
                        slope,amp,res=plfit(xvals,spec,[x0,x1])

                        if slope>0:
                            continue    

                        Line = 10**( slope*np.log10(xvals)+np.log10(amp))
                        chi_2 = ((np.log10(Line)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)

                        spectra_dict[LOS][sim].slope_bucket.slope_array[field][n0,n1] = slope
                        spectra_dict[LOS][sim].slope_bucket.chi_array[field][n0,n1] = chi_2
                        spectra_dict[LOS][sim].slope_bucket.slopext[field](nar([slope]))
                        slope_master(nar([slope]))
                        spectra_dict[LOS][sim].slope_bucket.chiext[field](nar([chi_2]))

if 0:
    #plot spectra
    for nf,field in enumerate(field_list):

        #field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
        if field in ['avg_cltt','avg_clee','avg_clbb']:
            ylim = [1e-11,1e6]
        else:
            ylim = [1e-10,0.5]
        #ylim = {'p':[1e-10,0.5], 'd':[1e-11,1e6]}

        for ns,sim in enumerate(sim_colors.simlist):
            if sim != 'half_half':
                continue
            print(field,sim)
            #for LOS in 'xyz':
            for LOS in 'y':
                if sim in ['avg_d','avg_v','avg_h'] and LOS in 'xy':
                    continue
                proj=rs.proj_dict[LOS][sim]
                x0_list = proj.lcent[2:10]
                x1_list = proj.lcent[11:51]
                Parr=np.zeros([x0_list.size,x1_list.size])
                Sarr=np.zeros([x0_list.size,x1_list.size])
                Chiarr=np.zeros([x0_list.size,x1_list.size])
                max_npoints=np.argmax(x1_list)-np.argmin(x0_list)+1
                fig,ax=plt.subplots(1,1)
                spec=spectra_dict[LOS][sim].spectra[field]

                ax.plot( proj.lcent, spec)

                if 1:
                    counter=0
                    for n0,x0 in enumerate(x0_list):
                        for n1,x1 in enumerate(x1_list):

                            ax2.clear()
                            ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                            xvals = proj.lcent[ok]
                                
                            spec=spectra_dict[LOS][sim].spectra[field][ok]
                            if 0:
                                #
                                # Linear Fit
                                #
                                slope,amp,res=plfit(xvals,spec,[x0,x1])
                            if 1:
                                #
                                # Quadratic Fit
                                #


                            if slope>0:
                                continue    

                            Line = 10**( slope*np.log10(xvals)+np.log10(amp))
                            chi_2 = ((np.log10(Line)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                            #ax.plot( xvals, Line, c=[0.5]*3 )
                            #ax2.plot(xvals, spec)
                            
                            #dt.axbonk(ax2,xscale='log',yscale='log',ylim=ylim,xlabel='k', xlim = [proj.lcent.min(), proj.lcent.max()])
                            #ax.plot( xvals, Line, c=[0.5]*3 )
                            #print(slope,do_ks_string( Line, spec), chi_2)
                            Kstat, crit, Pval = do_ks( Line, spec)
                            Parr[n0,n1]=Pval
                            Sarr[n0,n1]=slope
                            Chiarr[n0,n1]=chi_2


                            #x0,x1=proj.fitrange
                            #ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                            #xvals = proj.lcent[ok]
                            #slope = spectra_dict[LOS][sim].slopes[field]
                            #amp = spectra_dict[LOS][sim].amps[field]
                            #Line = 10**( slope*np.log10(xvals)+np.log10(amp))
                            #ax2.plot( xvals, Line, c='r',label="%0.2f"%slope)
                            #fig2.savefig('%s/OOO_%d.png'%(plotdir,counter))
                            #counter+=1


                #"Fit" line
                x0,x1=proj.fitrange
                ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                xvals = proj.lcent[ok]
                slope = spectra_dict[LOS][sim].slopes[field]
                amp = spectra_dict[LOS][sim].amps[field]
                Line = 10**( slope*np.log10(xvals)+np.log10(amp))
                ax.plot( xvals, Line, c='r',label="%0.2f"%slope)


                #<Fit> line

                #SA = spectra_dict[LOS][sim].slope_bucket.slope_array[field]
                weight = Parr
                slope = (Sarr*weight).sum()/weight.sum()

                do_it = not hasattr(spectra_dict[LOS][sim].slope_bucket, 'r_slope')
                if do_it:
                    spectra_dict[LOS][sim].slope_bucket.r_slope = {}
                spectra_dict[LOS][sim].slope_bucket.r_slope[field]=slope
                do_it = not hasattr(spectra_dict[LOS][sim].slope_bucket, 'p_array')
                if do_it:
                    spectra_dict[LOS][sim].slope_bucket.p_array = {}
                spectra_dict[LOS][sim].slope_bucket.p_array[field]=Parr
                do_it = not hasattr(spectra_dict[LOS][sim].slope_bucket, 's_array2')
                if do_it:
                    spectra_dict[LOS][sim].slope_bucket.s_array2 = {}
                spectra_dict[LOS][sim].slope_bucket.s_array2[field]=Sarr

                do_it = not hasattr(spectra_dict[LOS][sim].slope_bucket, 'chi_array2')
                if do_it:
                    spectra_dict[LOS][sim].slope_bucket.chi_array2 = {}
                spectra_dict[LOS][sim].slope_bucket.chi_array2[field]=Chiarr

                #inv_chi=1./spectra_dict[LOS][sim].slope_bucket.chi_array[field]
                #tots=inv_chi.sum()
                #slope = (SA*inv_chi).sum()/tots
                amp = spectra_dict[LOS][sim].slope_bucket.amp16[field]
                xvals = proj.lcent
                Line = 10**( slope*np.log10(xvals))
                Line = Line*amp/Line[16]

                ax.plot( xvals, Line, c='g',label="%0.2f"%slope)
                ax.legend(loc=0)

                dt.axbonk(ax,xscale='log',yscale='log', ylim=ylim)
                fig.savefig('%s/Specs_%s_%s_%s.png'%(plotdir,field,LOS,sim))
                plt.close(fig)


if 1:
    plt.close('all')
    got_one=False
    for LOS in 'y':#'xyz':
        for nf,field in enumerate(field_list):
            if got_one:
                break
            fig,ax=plt.subplots(1,2)
            for ns,sim in enumerate(sim_colors.simlist):
                do_it = not hasattr(spectra_dict[LOS][sim].slope_bucket, 'r_slope')
                if do_it:
                    continue
                chi_arr=spectra_dict[LOS][sim].slope_bucket.chi_array2[field].flatten()
                ok = chi_arr != 0
                chi_arr = chi_arr[ok]
                log_chi = np.log10(chi_arr)
                chi_map = dt.rainbow_map( chi_arr.max() - chi_arr.min())
                got_one=True
                print("XXXXXXXX")
                Sarr=spectra_dict[LOS][sim].slope_bucket.s_array2[field].flatten()[ok]
                Parr=spectra_dict[LOS][sim].slope_bucket.p_array[field].flatten()[ok]
                #print(Sarr)
                #print(Parr)
                c=chi_map( chi_arr.flatten() - chi_arr.min())
                ax[0].scatter( Parr.flatten(), Sarr.flatten(), c=c)
                norm = mpl.colors.Normalize(vmin=0,vmax=1)
                line=ax[1].scatter( chi_arr, Sarr, c=Parr,norm=norm,s=1)
                ax[1].set_xscale('log')
                fig.colorbar(line)
                #proj=rs.proj_dict[LOS][sim]
                #S1 = spectra_dict[LOS][sim].slopes[field]
                #Sarr = spectra_dict[LOS][sim].slope_bucket.slope_array[field]
                #S2 = Sarr.mean()
                #V2 = Sarr.std()
                #ax.errorbar(S1, S2, yerr=V2,marker=sim_colors.marker[sim],color=sim_colors.color[sim])
            #ax.plot([-4.5,-2.5],[-4.5,-2.5],c=[0.5]*4)
            outname='%s/SlopeVsKSP_%s_%s.png'%(plotdir,field,LOS)
            fig.savefig(outname)
            print(outname)
if 0:
    plt.close('all')
    for LOS in 'xyz':
        for nf,field in enumerate(field_list):
            fig,ax=plt.subplots(1,1)
            for ns,sim in enumerate(sim_colors.simlist):
                print(field,sim)
                proj=rs.proj_dict[LOS][sim]
                S1 = spectra_dict[LOS][sim].slopes[field]
                Sarr = spectra_dict[LOS][sim].slope_bucket.slope_array[field]
                S2 = Sarr.mean()
                V2 = Sarr.std()
                ax.errorbar(S1, S2, yerr=V2,marker=sim_colors.marker[sim],color=sim_colors.color[sim])
            ax.plot([-4.5,-2.5],[-4.5,-2.5],c=[0.5]*4)
            fig.savefig('%s/SlopeVsSlope_%s_%s.png'%(plotdir,field,LOS))

if 0:
    plt.close('all')
    for LOS in 'xyz':
        for nf,field in enumerate(field_list):
            for ns,sim in enumerate(sim_colors.simlist):
                if field in ['avg_d','avg_v','avg_h'] and LOS in 'xz':
                    continue
                do_it = not hasattr(spectra_dict[LOS][sim].slope_bucket, 'r_slope')
                if do_it:
                    continue
                fig,ax=plt.subplots(1,1)
                print(field,sim)
                proj=rs.proj_dict[LOS][sim]
                S1 = spectra_dict[LOS][sim].slopes[field]
                Sarr = spectra_dict[LOS][sim].slope_bucket.slope_array[field]
                Chiarr = spectra_dict[LOS][sim].slope_bucket.chi_array[field]
                InvChi=1./Chiarr
                S2 = Sarr.mean()
                V2 = Sarr.std()
                hist,bins,obj=ax.hist( Sarr.flatten(), histtype='step')
                ax.scatter( S1, 1.2*hist.max(),marker=sim_colors.marker[sim],color=sim_colors.color[sim])

                ax.scatter( (Sarr*InvChi).sum()/InvChi.sum() , 1.3*hist.max(),marker=sim_colors.marker[sim],color='k')


                ax.errorbar(S2, hist.max(), xerr=V2,marker=sim_colors.marker[sim],color=sim_colors.color[sim])
                
                ax.errorbar(spectra_dict[LOS][sim].slope_bucket.r_slope[field], hist.max(), c='r', marker='s')


                lcent=proj.lcent
                spec=spectra_dict[LOS][sim].spectra[field]
                slope,amp,res=plfit(lcent,spec,fitrange)

                ax.scatter( slope, 2*hist.max(),c='g',marker='*')
                #dt.axbonk(ax,xlim=slope_master.minmax)






                fig.savefig('%s/HistVsSlope_%s_%s_%s.png'%(plotdir,field,LOS,sim))
                plt.close(fig)

if 0:
    plt.close('all')

    counter=0
    slope_old=[]
    slope_new=[]    
    f3,ax3=plt.subplots(1,1)
    #ax2=ax2.flatten()

    #for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb'])
    chiext=dt.extents()
    slopext=dt.extents()
    for nf,field in enumerate(['avg_clbb']):
        #for aaa in ax2:
        #    aaa.clear()
        ax3.clear()
        for ns,sim in enumerate(sim_colors.simlist):
            f2,ax2=plt.subplots(2,2)
            ax2=ax2.flatten()
            counter += 1
            print(counter)
            #if counter not in [2]:
            #    continue
            #if sim != '1_2':
            #    continue

            proj=rs.proj_dict['x'][sim]
            label="%s %0.1f"%(sim, spectra_dict['x'][sim].slopes[field])

            x0_list = proj.lcent[2:10:2]
            x1_list = proj.lcent[11:41:4]
            max_npoints=np.argmax(x1_list)-np.argmin(x0_list)+1
            rm = dt.rainbow_map(max_npoints)
            rm0 = dt.rainbow_map(len(x0_list))
            rm = dt.rainbow_map(len(x1_list))

            fig4,ax4=plt.subplots(1,1)
            ntot=0
            spectra_dict[LOS][sim].spext=dt.extents()
            spectra_dict[LOS][sim].speyt=dt.extents()
            xlim = {'1_half':[3.07e-02, 4.85e-01], '1_2':[3.07e-02, 4.85e-01]}.get(sim,None)

            ylim = {'1_half':[1.77e-05, 5.93e-01], '1_2':[4.73e-03, 1.03e+02]}.get(sim,None)
            slopelim=  [-4.75e+00, -1.81e+00]
            chilim = [1.24e-04, 2.05e+03]

            spectra_dict[LOS][sim].slope_array= np.zeros([x0_list.size, x1_list.size]) +10
            spectra_dict[LOS][sim].chi_array= np.zeros([x0_list.size, x1_list.size]) +10
            for n0,x0 in enumerate(x0_list):
                for n1,x1 in enumerate(x1_list):

                    ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                    xvals = proj.lcent[ok]
                        
                    #axlist[nf].plot(proj.lcent[ok],   spectra_dict['x'][sim].spectra[field][ok], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)

                    spec=spectra_dict['x'][sim].spectra[field][ok]
                    slope,amp,res=plfit(xvals,spec,[x0,x1])

                    Line = 10**( slope*np.log10(xvals)+np.log10(amp))
                    #axlist[nf].plot( xvals, Line,c='k',label='F1')
                    spectra_dict[LOS][sim].slope_array[n0,n1] = slope
                    if slope>0:
                        continue    

                    #KS_stat, KS_crit, KS_p=do_ks( np.log10(Line), np.log10(spec))
                    ##ax2.scatter( slope, ( (Line-spec)**2/spec).sum(), c=[rm( ok.sum())])
                    #chi_2 = ((Line-spec)**2/spec).sum()/( ok.sum()-2)
                    chi_2 = ((np.log10(Line)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)

                    spectra_dict[LOS][sim].chi_array[n0,n1] = chi_2
                    spectra_dict[LOS][sim].slopext(nar([slope]))
                    spectra_dict[LOS][sim].chiext(nar([chi_2]))
                    #continue
                    ax2[0].scatter( slope, chi_2, c=[rm( ok.sum())])
                    #ax2[0].scatter( slope, chi_2, c=[rm(n0)])
                    ax2[2].scatter( ok.sum(), chi_2,c='k')
                    #ax2[2].scatter( ok.sum(), chi_2/ok.sum(), c='r')
                    #ax2[3].scatter( ok.sum(), slope, c=[rm(n1)])
                    ax2[3].scatter( ok.sum(), slope, c=[rm0(n0)])

                    ax4.clear()
                    ax4.plot( xvals, Line, c='g')
                    ax4.plot( xvals, spec, c='k')
                    spext(xvals);speyt(Line); speyt(spec)
                    dt.axbonk(ax4,xlim=xlim, ylim=ylim,xscale='log',yscale='log')
                    fig4.savefig('%s/OOF_%s_%d.png'%(plotdir,sim,ntot))
                    ntot+=1 


            if 0:
                n = spec.spectra[field].size
                alphahat = 1 + ok.sum()/( np.log(spec.spectra[field][ok]/fitrange[0]).sum())
                print("SLOOOO", spec.slopes[field], alphahat)
                dk = proj.lbins[1:]-proj.lbins[:-1]
                A = np.sum( spec.spectra[field][ok]*dk[ok])*(1-alphahat)*fitrange[0]**alphahat/(fitrange[1]**(1-alphahat)-fitrange[0]**(1-alphahat))

                Line2 = 10**( alphahat*np.log10(proj.lcent))
                Line2 *= spec.spectra[field][ok][5]/Line2[ok][5]
                axlist[nf].plot( proj.lcent[ok], Line2[ok],label='F2',c='r')
            #bins = np.mgrid[ slope_array.min():slope_array.max():16j]
            bins = np.mgrid[ -4.5:0:16j] +np.random.random()*0.01   
            bins = np.mgrid[slope_array.min(): slope_array.min()+4.5:0.3]
            print(bins)
            #if slope_array.min()<bins.min():
            #    raise
            #if slope_array.max()>bins.max():
            #    raise

            #ax2[1].hist( slope_array[slope_array < 0].flatten())
            ax3.clear()
            ax3.hist( slope_array[slope_array < 0].flatten(), histtype='step', linestyle=sim_colors.line_list[ns],color=sim_colors.color_list[ns],bins=bins)
            dt.axbonk(ax3,xlabel=field,ylabel='N',xlim=[-4.5,0],ylim=[0,200])
            f3.savefig('%s/refit_SlopeHist_%s_%s.png'%(plotdir,field,sim))

            ax2[1].hist( slope_array[slope_array < 0].flatten(), histtype='step', linestyle=sim_colors.line_list[ns],color=sim_colors.color_list[ns],bins=bins)
            dt.axbonk(ax2[1],xlabel=field,ylabel='N',xlim=[-4.5,0],ylim=[0,200])
            dt.axbonk(ax2[0],xlabel='slope',ylabel='Chi',yscale='log', ylim=chilim,xlim=slopelim)
            dt.axbonk(ax2[2],xlabel='Ndof',ylabel='Chi',yscale='log', ylim=chilim)
            dt.axbonk(ax2[3],xlabel='Ndof',ylabel='Slope', ylim=slopelim)


            #proj=rs.proj_dict['y'][sim]
            #label="%s %0.1f"%(sim, spectra_dict['y'][sim].slopes[field])
            #axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim], label=label)
            #axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            f2.savefig('%s/refit_SlopeThings_%s_%s.png'%(plotdir,field,sim))
            plt.close(f2)

