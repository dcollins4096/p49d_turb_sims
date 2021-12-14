
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

from scipy.optimize import curve_fit
def cube(x,a,b,c,d):
    return a + b*x + c*x**2 +d*x**3
def cube_prime(x,a,b,c,d):
    return  b + 2*c*x +3*d*x**2
def cube_prime_prime(x,a,b,c,d):
    return   2*c +6*d*x
def quad(x,a,b,c):
    return a + b*x + c*x**2
def line(x,a,b):
    return a + b*x

def get_zero_substrings(arr,thresh=1e-16):
    starts = [0]
    lengths = [0]
    for i, v in enumerate(arr):
        if np.abs(v) < thresh:
            lengths[-1]+=1
        else:
            starts.append(i)
            lengths.append(0)
    return starts, lengths


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
        self.quad_fits={}
        self.lin_fits={}

field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
simlist = ['3_half'] # ['half_half']#,'2_2']
field_list = ['avg_clbb']
#   if 0:
#       ax_rho_p.clear()
#       from scipy.optimize import curve_fit
#       def density_ac(r, A, r0, gamma):
#           return A*(1+(r/r0)**gamma)
#       popt, pcov = curve_fit(density_ac, the_x, AC, p0=[1,1,-2])
#       po2 = AC[0]/2, 1, -1
#       the_x = np.linspace(0,0.25,100)
#       ax_rho_p.plot( the_x, density_ac(the_x, po2[0],po2[1],po2[2]),c='k')

def compensate_plot(ax,alpha,x,y,**kwargs):
    ax.plot( x, x**alpha*y, **kwargs)
if 1:
    #
    # Development code.  Now part of read_stuff
    #
    #fill objects with slope array
    slope_master = dt.extents()
    spec_master = dt.extents()
    counter=0
    for ns,sim in enumerate(sim_colors.simlist):
        spectra_dict[LOS][sim].slope_bucket = slope_bucket()

    fig_all,ax_all=plt.subplots(1,2, figsize=(8,12))
    for nf,field in enumerate(field_list):
        ax_all[0].clear()
        ax_all[1].clear()
        bbb=[]
        ccc=[]
        sss=[]
        sssm=[]
        MS=[]
        MA=[]
        aaa=[]
        LOS_LIST = 'y'
        if field in ['avg_d','avg_v','avg_h']:
            LOS_LIST='x'
        for ns,sim in enumerate(simlist):
            
            for LOS in LOS_LIST:
                counter+=1
                #print(field,sim, LOS, counter)
                #continue

                #if counter not in [1, 63]:
                #    continue
                proj=rs.proj_dict[LOS][sim]
                specf=spectra_dict[LOS][sim].spectra[field]
                fit_range=proj.determine_fit_range()
                x0=fit_range[0]
                x1=fit_range[1]

                ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                xvals = proj.lcent[ok]
                spec=spectra_dict[LOS][sim].spectra[field][ok]
                popt_fixed, pcov_fixed= curve_fit( line, np.log10(xvals), np.log10(spec))

                AlphaF = popt_fixed[1]

                fig,axes=plt.subplots(2,2,figsize=(12,4))
                axL = axes.flatten()

               # axL[1].plot([-12,3],[-12,3],c=[0.5]*3)
                axL[1].plot([-3,0],[-3,0],c=[0.5]*3)
                axL[3].plot([1e-6,200],[1e-6,200],c='k')
                axL[3].plot([1e-6,200],[0.1*1e-6,0.1*200],c='k')
                axL[3].plot([1e-6,200],[10*1e-6,10*200],c='k')
                #axL[3].set_yscale('log')
                #dt.axbonk(axL[2],xlabel=r'$\alpha_3$', ylabel=r'$\chi^2_3$', xlim=[-3,0], ylim=[1e-3,1e2])
                #axL[2].set_yscale('log')


                x0_list = [x0]
                x1_list = [x1]
                x0_list = proj.lcent[2:10:2]
                x1_list = proj.lcent[12:51:8]

                counter2=0
                chi_l_list = []
                slope_list = []
                amp_list = []

                slope_list_3=[]
                amp_list_3=[]
                k_list_3=[]
                ok_list=[]
                color_3=[]
                a_list=[]
                b_list=[]
                c_list=[]
                d_list=[]

                norm_3 = mpl.colors.Normalize(vmin=min(x0_list), vmax=max(x0_list))
                norm_3 = mpl.colors.LogNorm( vmin=1e-3, vmax=1)
                norm_dof = mpl.colors.Normalize(vmin=1,vmax=50)
                for n0,x0 in enumerate(x0_list):
                    for n1,x1 in enumerate(x1_list):
                        color_3.append(n0)
                        ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                        xvals = proj.lcent[ok]
                        spec=spectra_dict[LOS][sim].spectra[field][ok]
                        popt_3, pcov_3= curve_fit( cube, np.log10(xvals), np.log10(spec))
                        popt_1, pcov_1= curve_fit( line, np.log10(xvals), np.log10(spec))


                        this_c=10**cube(np.log10(xvals), popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                        this_l=10**line(np.log10(xvals), popt_1[0],popt_1[1])

                        chi_c = ((np.log10(this_c)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                        chi_l = ((np.log10(this_l)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)

                        chi_l_list.append(chi_l)
                        slope_list.append(popt_1[1])
                        amp_list.append(popt_1[0])
                        ok_list.append(ok)

                        #Kstat1, crit1, Pval1 = do_ks( np.log10(this_l), np.log10(spec))
                        #Kstat2, crit2, Pval2 = do_ks( np.log10(this_q), np.log10(spec))

                        a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                        a_list.append(a3);b_list.append(b3); c_list.append(c3); d_list.append(d3)
                        slope = b3 - c3**2/(3*d3)
                        value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                        k_flat = -c3/(3*d3)
                        slope_list_3.append(slope)
                        amp_list_3.append(value)
                        k_list_3.append(k_flat)

                        #line_from_cube = 10**( AlphaF*np.log10(xvals-k_flat)+value)
                        line_from_cube =10**line(np.log10(xvals), popt_1[0],popt_1[1])
                        line_from_cube = 10**( slope*(np.log10(xvals)-k_flat) + value)


                        if 0:
                            import ks_b
                            F1=np.log10(this_l)
                            
                            F3=np.log10(spec)
                            #d, prob, ne = ks_b.kstest( F1, F2)
                            axL[3].clear()
                            axL[3].plot( sorted(F1), np.arange(0,1, 1./F1.size),c='g')
                            axL[3].plot( sorted(F2), np.arange(0,1, 1./F2.size),c='r')
                            axL[3].plot( sorted(F3), np.arange(0,1, 1./F3.size),c='k')

                        axL[0].clear()
                        axL[0].set_title(r'$\alpha_1 = %0.2f \alpha_3 = %0.2f \chi^2_1 = %0.1e \chi^2_3 = %0.1e$'%(popt_1[1], slope, chi_l, chi_c))
                        #axL[0].plot( proj.lcent, specf, c=[0.5]*4)
                        compensate_plot(axL[0], -AlphaF, proj.lcent, specf, c=[0.5]*4)
                        compensate_plot(axL[0], -AlphaF, xvals, this_c, c='b')
                        #compensate_plot(axL[0], -AlphaF, xvals, this_q, c='r')
                        compensate_plot(axL[0], -AlphaF,  xvals, this_l,c='g')
                        compensate_plot(axL[0], -AlphaF,  xvals, line_from_cube,c='r')
                        if field in ['avg_cltt','avg_clee','avg_clbb']:
                            #ylim = [1e-11,1e6] #whole range
                            ylim = [1e-3, 1e-2]
                        else:
                            ylim = [1e-10,0.5]
                            #ylim = [1e-4, 1e-1]
                            ylim = [1e-5, 1e-3]
                        dt.axbonk(axL[0],xscale='log',yscale='log', ylim=ylim)

                        #print("Fix %0.2e Cube %0.2e"%(popt_1[1], slope))
                        #axL[1].scatter( popt_1[1],slope, c=x0, norm=norm_3)
                        axL[1].scatter( popt_1[1],slope, c=chi_c, norm=norm_3)
                        dt.axbonk(axL[1], xlim=[-3,0],ylim=[-3,0], xlabel=r'$\alpha_1$', ylabel=r'$\alpha_3$')
                        #dt.axbonk(axL[1],xlim=[-12,3],ylim=[-12,3], xlabel=r'$\alpha_1$',ylabel=r'$\alpha_2$')
                        if 0:
                            axL[1].scatter( ok.sum(), popt_1[0],c='r',marker='o')
                            axL[1].scatter( ok.sum(), popt_1[1],c='r',marker='s')
                            axL[1].scatter( ok.sum(), popt_2[0],c='g',marker='o')
                            axL[1].scatter( ok.sum(), popt_2[1],c='g',marker='s')
                            axL[1].scatter( ok.sum(), popt_2[2],c='g',marker='v')

                            axL[1].scatter( ok.sum(), popt_3[0],c='b',marker='o')
                            axL[1].scatter( ok.sum(), popt_3[1],c='b',marker='s')
                            axL[1].scatter( ok.sum(), popt_3[2],c='b',marker='v')
                            axL[1].scatter( ok.sum(), popt_3[3],c='b',marker='^')


                        if 0:
                            this_cprime=10**cube_prime(np.log10(xvals), popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                            this_cprimeprime=10**cube_prime_prime(np.log10(xvals), popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                            #makes useful plots.
                            #axL[2].clear()
                            axL[2].plot( xvals, this_cprime,c='r')
                            axL[2].plot( xvals, this_cprimeprime,c='g')
                            #axL[2].set_xscale('log'); axL[2].set_xlim( axL[0].get_xlim())
                            #axL[2].set_yscale('log')

                        if 1:
                            MYX = np.log10(xvals)
                            this_c=cube(MYX, popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                            this_cprime=cube_prime(MYX, popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                            this_cprimeprime=cube_prime_prime(MYX, popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                            a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                            k_flat = -2*c3/(6*d3)
                            slope = b3 - c3**2/(3*d3)
                            value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                            slope_peak= -2*c3/(6*d3)
                            #index = np.argmin( np.abs(MYX -slope_peak))
                            #slope1= cube_prime( slope_peak, popt_3[0],popt_3[1],popt_3[2],popt_3[3])
                            slope = b3 - c3**2/(3*d3)
                            value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)

                            #value = cube( slope_peak, a3,b3,c3,d3)
                            #print(value, ampl)

                            #axL[2].scatter( slope, chi_c, c=ok.sum(),norm=norm_dof)
                            #print(ok.sum())
                            #axL[2].plot( MYX, this_c, c=[0.5]*3)
                            #axL[2].scatter( slope_peak, value)
                            #axL[2].scatter( MYX, (MYX-slope_peak)*slope+value)
                            #print("NEW",slope)
                            #new_line = 10**( (xvals-slope_peak)*slope+value)
                            #axL[2].plot( MYX, this_cprime,c='k')
                            #axL[2].plot( MYX, a3+b3*MYX+c3*MYX**2+d3*MYX**3)
                            #axL[2].plot( MYX, cube_prime(MYX, a3,b3,c3,d3))
                            #axL[2].plot( MYX, b3+2*c3*MYX+3*d3*MYX*2)
                            #axL[2].plot( xvals, this_c, c=[0.5]*3)
                            #axL[2].plot( xvals, new_line)
                            #axL[2].set_xscale('log'); axL[2].set_xlim( axL[0].get_xlim())
                            #axL[2].set_yscale('log')



                        #print(chi_l, chi_q)
                        axL[3].scatter( chi_l, chi_c)
                        dt.axbonk(axL[3],xscale='log',yscale='log', xlim=[1e-6,200],ylim=[1e-6,200], xlabel=r'$\chi^2_1$',ylabel=r'$\chi^3_2$')


                        #axL[3].scatter( Pval1, Pval2)
                        #print(Pval1, Pval2)
                        #axL[3].scatter(popt_2[2]/popt_1[1], chi_l)
                        #axL[3].scatter(popt_1[1], popt_2[2]/popt_1[1],norm=norm,c=chi_l)
                        axL[3].scatter( popt_1[1], chi_l )

                        counter2+=1
                        fig.savefig('%s/quad_%s_%s_%s_%d.png'%(plotdir,field,sim,LOS, counter2))
                        


                plt.close('all')
                import three_way_bean
                plt.clf()
                fig4, ax_alpha, ax_top, ax_right = three_way_bean.three_way_bean()

                bins = np.mgrid[-5:0:16j]
                norm = mpl.colors.Normalize(vmin=min(color_3), vmax=max(color_3))
                ax_alpha.scatter( slope_list, slope_list_3, c=color_3,norm=norm)
                ax_alpha.plot( [-5,0],[-5,0],c=[0.5]*3)
                hist, bins, visthings=ax_right.hist( slope_list_3, histtype='step', orientation='horizontal',bins=bins)
                ax_top.hist( slope_list, histtype='step',bins=bins)
                dt.axbonk(ax_alpha, xlabel=r'$\alpha_1$',ylabel=r'$\alpha_3$', xlim=[-5,0],ylim=[-5,0])

                ax_right.set_ylim( ax_alpha.get_ylim())
                ax_top.set_xlim( ax_alpha.get_xlim())

                most_probable_bin = np.argmax(hist)
                bL = bins[most_probable_bin]
                bR = bins[most_probable_bin+1]
                SL = nar(slope_list_3)
                most_probable =  (SL >= bL)*(SL <= bR)
                AlphaC = SL[most_probable].mean()
                ax_alpha.axvline(AlphaC, c=[0.5]*3)
                ax_alpha.axhline(AlphaC, c=[0.5]*3)
                AL = nar(amp_list_3)
                Amp = AL[ most_probable ].mean()
                KL = nar(k_list_3)
                Kfit = KL[ most_probable ].mean()

                fig4.savefig('%s/alpha3_%s_%s_%s_%d.png'%(plotdir,field,sim,LOS, counter2))
                plt.close(fig4)

                fig2,ax2=plt.subplots(1,1, figsize=(12,8))
                x0=fit_range[0]; x1=fit_range[1]
                ok = (proj.lcent>=x0)*(proj.lcent<=x1)
                xvals = proj.lcent[ok]
                compensate_plot(ax2,-AlphaC, proj.lcent, specf, c=[0.5]*4)


                #
                # Fixed
                #

                this_l=10**line(np.log10(xvals), popt_1[0],AlphaF)
                kfit = proj.lcent[16]
                index = np.argmin( np.abs( xvals-kfit))
                this_l = this_l*specf[16]/this_l[index]
                chi2 = np.sum( (this_l-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, this_l, c='g',label=r'$%0.2e \pm %0.2e$'%(AlphaF,chi2))
                


                #
                # Mean
                #

                line_from_cube = 10**( AlphaC*(np.log10(xvals)-Kfit) + Amp)
                chi2 = np.sum( (line_from_cube-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, line_from_cube, c='r',label=r'$%0.2e \pm %0.2e$'%(AlphaC,chi2))

                #
                # cube
                #

                #most_probable = slice(3,4)
                mean_a = nar(a_list)[most_probable].mean()
                mean_b = nar(b_list)[most_probable].mean()
                mean_c = nar(c_list)[most_probable].mean()
                mean_d = nar(d_list)[most_probable].mean()
                sss.append( AlphaC)
                sssm.append( mean_b - mean_c**2/mean_d)
                a_list=nar(a_list)
                b_list=nar(b_list)
                c_list=nar(c_list)
                d_list=nar(d_list)


                #index = np.argmin( np.abs( xvals-Kfit))
                #cube_from_cube = 10**( AlphaC*(np.log10(xvals)-Kfit) + Amp)
                #cube_from_mean = 10**cube( np.log10(xvals), popt_3[0], popt_3[1], popt_3[2], popt_3[3])
                cube_from_mean = 10**cube( np.log10(xvals), mean_a, mean_b, mean_c, mean_d)
                #ax2.scatter( xvals[index],xvals[index]**(-AlphaC)*cube_from_mean[index])
                #ax2.scatter( 10**Kfit, cube_from_mean.min())
                chi2 = np.sum( (cube_from_mean-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, cube_from_mean, c='b',label=r'$%0.2e \pm %0.2e$'%(AlphaC,chi2))

#                if field in ['avg_cltt','avg_clee','avg_clbb']:
#                    ylim = [1e-11,1e6] #whole range
#                    #ylim = [1e-3, 1e-2]
#                else:
#                    ylim = [1e-10,0.5]
                ax2.legend(loc=0)
                dt.axbonk( ax2, xscale='log',yscale='log')#, ylim=ylim)
                #print('%s/Net_%s_%s_%s.png'%(plotdir,field,sim,LOS))
                #fig2.savefig('%s/Net_%s_%s_%s.png'%(plotdir,field,sim,LOS))

                ax_all[0].scatter( mean_a, mean_c**2/(3*mean_d), c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( AlphaC, -mean_c**2/(3*mean_d), c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( c_list/b_list,d_list/b_list, c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( c_list[most_probable]/b_list[most_probable],d_list[most_probable]/b_list[most_probable], c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                ax_all[1].scatter( c_list[most_probable],d_list[most_probable], c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( mean_c/mean_b,mean_d/mean_b, c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                bbb.append(mean_a); ccc.append(mean_c**2/(3*mean_d))
                SL2 = nar(c_list)**2/(3*nar(d_list))
                #ax_all[0].scatter( b_list[most_probable], SL2[most_probable],c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[0].scatter( b_list, SL2,c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                print(sim.split("_"), quan3[sim]['msavg'])
                #ax_all[1].scatter( quan3[sim]['msavg'], mean_a, c=sim_colors.color[sim],marker=sim_colors.marker[sim])


                MS.append( quan3[sim]['msavg'])
                MA.append( quan3[sim]['maavg'])
                aaa.append(mean_a)

        def fitFunc(x, a, b, c, d):
            return a+ b*x[0] + c*x[1] + d*x[0]*x[1]
        def do_stuff(xdata, ydata, p0):
            fitParams, fitCovariances = curve_fit(fitFunc, xdata, ydata, p0)   
            print(' fit coefficients:\n', fitParams)
            return fitParams, fitCovariances
        def do_morestuff(aax, myx, myy, myval,c='k'):
            x=np.row_stack([myx.flatten(),myy.flatten()])
            y=myval.flatten()
            p0 = [0,1,1,0]
            fitParams, fitCovariances= do_stuff(x,y, p0)
            f = fitParams
            amps = f[0] + f[1]*myx + f[2]*myy+f[3]*myx*myy
            chi2 = ( (myval-amps)**2/amps).sum()
            print("chi", c,chi2)
            aax.scatter(myx, amps, c=c,s=0.2)
        #do_morestuff( ax_all[1], nar(MS),nar( MA), nar(aaa))

        minmin=min([min(bbb),min(ccc)])
        maxmax=max([max(bbb),max(ccc)])
        bbb=nar(bbb)
        pfit = np.polyfit( bbb,ccc,1)

        ax_all[0].plot( bbb, bbb*pfit[0]+pfit[1])
        print(bbb)
        #ax_all[0].plot([minmin,maxmax],[minmin,maxmax],c=[0.5]*3)
        dt.axbonk(ax_all[0],xlabel='b',ylabel=r'$c^2/(3 d)$')
        fig_all.savefig('%s/Linearity_%s.png'%(plotdir,field))
        print("FARTS %s"%field)

                #ax.plot( proj.lcent, specf, c=[0.5]*4)
                #ax.plot( xvals, 10**quad(np.log10(xvals), popt_2[0],popt_2[1],popt_2[2]), c='r')
                #ax.plot( xvals, 10**line(np.log10(xvals), popt_1[0],popt_1[1]),c='g')
                #dt.axbonk(ax,xscale='log',yscale='log')
                #fig.savefig('%s/quad_%s_%s_%s.png'%(plotdir,field,sim,LOS))

