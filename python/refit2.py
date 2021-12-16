field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
simlist = ['half_2', 'half_half']#['half_half']#['1_1']#['2_2'] #['half_half'] # ['half_half']#,'2_2']
field_list = ['avg_cltt']
do_all_subplots=False

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
def minmax(alpha,x,y):
    v = x**alpha*y
    return nar([min(v),max(v)])
                        

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
from scipy.ndimage import gaussian_filter1d
def compensate_plot(ax,alpha,x,y,**kwargs):
    ax.plot( x, x**alpha*y, **kwargs)
if 1:
    #
    # Development code.  
    #
    #fill objects with slope array
    for ns,sim in enumerate(sim_colors.simlist):
        spectra_dict[LOS][sim].slope_bucket = slope_bucket()

    fig_all,ax_all=plt.subplots(1,1, figsize=(8,12))
    for nf,field in enumerate(field_list):
        #ax_all[0].clear()
        #ax_all[1].clear()
        bbb=[]
        ccc=[]
        MS=[]
        MA=[]
        LOS_LIST = 'y'
        if field in ['avg_d','avg_v','avg_h']:
            LOS_LIST='x'

        for ns,sim in enumerate(simlist):
            
            for LOS in LOS_LIST:
                print("DO", field, sim, LOS)
                proj=rs.proj_dict[LOS][sim]
                lcentf = proj.lcent
                specf=spectra_dict[LOS][sim].spectra[field]
                if 1:
                    fit_range=proj.determine_fit_range()
                    x0_list = lcentf[2:10:2]
                    x1_list = lcentf[14:51:8]
                if 0:
                    #fit_range=proj.determine_fit_range()
                    fit_range = proj.lcent[2], proj.lcent[23] #two less than the fit range for dx1.  queb3 line 324
                    specf = gaussian_filter1d( specf[2:], 1, mode='nearest')
                    lcentf = lcentf[2:]
                    x0_list = lcentf[0:10:2]
                    x1_list = lcentf[12:51:8]
                x0=fit_range[0]
                x1=fit_range[1]

                #x0_list = [x0]
                #x1_list = [x1]
                #x0_list = [2,3]
                #x1_list = [45,51]

                # fixed fit
                ok = (lcentf>=x0)*(lcentf<=x1)
                xvals = lcentf[ok]
                spec=specf[ok]
                popt_fixed, pcov_fixed= curve_fit( line, np.log10(xvals), np.log10(spec))
                AlphaF = popt_fixed[1]

                fig,axes=plt.subplots(2,3,figsize=(12,4))
                axL = axes.flatten()

                axL[2].plot([1e-6,200],[1e-6,200],c='k')
                axL[2].plot([1e-6,200],[0.1*1e-6,0.1*200],c='k')
                axL[2].plot([1e-6,200],[10*1e-6,10*200],c='k')



                counter2=0
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

                dof_list=[]
                verts = dt.extents()
                c_to_b=dt.extents()
                d_to_b=dt.extents()
                a3_ext=dt.extents()
                a1_ext=dt.extents()
                chi_c_ext=dt.extents()
                chi_l_ext=dt.extents()
                a_ext=dt.extents()
                b_ext=dt.extents()

                #
                # Loop for finding extents and filling arrays
                #
                for n0,x0 in enumerate(x0_list):
                    for n1,x1 in enumerate(x1_list):
                        ok = (lcentf>=x0)*(lcentf<=x1)
                        xvals = lcentf[ok]
                        spec=specf[ok]
                        popt_3, pcov_3= curve_fit( cube, np.log10(xvals), np.log10(spec))
                        popt_1, pcov_1= curve_fit( line, np.log10(xvals), np.log10(spec))

                        a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                        slope = b3 - c3**2/(3*d3)
                        value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                        k_flat = -c3/(3*d3)
                        if 10**k_flat <= x0 or 10**k_flat >= x1:
                            continue
                        this_c=10**cube(np.log10(xvals), popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                        this_l=10**line(np.log10(xvals), popt_1[0],popt_1[1])

                        chi_c = ((np.log10(this_c)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                        chi_l = ((np.log10(this_l)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                        chi_c_ext(chi_c)
                        chi_l_ext(chi_l)
                        c_to_b(c3/b3)
                        d_to_b(d3/b3)
                        a_ext(a3)
                        b_ext(b3)
                        a1_ext(popt_1[1])
                        a3_ext(slope)
                        slope_list.append(popt_1[1])
                        amp_list.append(popt_1[0])
                        ok_list.append(ok)
                        dof_list.append(ok.sum())
                        slope_list_3.append(slope)
                        amp_list_3.append(value)
                        k_list_3.append(k_flat)
                        color_3.append(n0)
                        a_list.append(a3);b_list.append(b3); c_list.append(c3); d_list.append(d3)
                #
                # Do a bunch of plots
                #
                for n0,x0 in enumerate(x0_list):
                    for n1,x1 in enumerate(x1_list):
                        if not do_all_subplots:
                            continue
                        ok = (lcentf>=x0)*(lcentf<=x1)
                        xvals = lcentf[ok]
                        spec=specf[ok]
                        popt_3, pcov_3= curve_fit( cube, np.log10(xvals), np.log10(spec))
                        popt_1, pcov_1= curve_fit( line, np.log10(xvals), np.log10(spec))

                        a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                        slope = b3 - c3**2/(3*d3)
                        value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                        k_flat = -c3/(3*d3)
                        if 10**k_flat <= x0 or 10**k_flat >= x1:
                            continue

                        this_c=10**cube(np.log10(xvals), popt_3[0],popt_3[1],popt_3[2], popt_3[3])
                        this_l=10**line(np.log10(xvals), popt_1[0],popt_1[1])

                        chi_c = ((np.log10(this_c)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)
                        chi_l = ((np.log10(this_l)-np.log10(spec))**2/spec).sum()/( ok.sum()-2)

                        line_from_cube =10**line(np.log10(xvals), popt_1[0],popt_1[1])
                        line_from_cube = 10**( slope*(np.log10(xvals)-k_flat) + value)

                        axL[0].clear()
                        axL[0].set_title(r'$\alpha_1 = %0.2f \alpha_3 = %0.2f \chi^2_1 = %0.1e \chi^2_3 = %0.1e$'%(popt_1[1], slope, chi_l, chi_c))

                        compensate_plot(axL[0], -AlphaF, lcentf, specf, c=[0.5]*4)
                        compensate_plot(axL[0], -AlphaF, xvals, this_c, c='b')
                        compensate_plot(axL[0], -AlphaF,  xvals, this_l,c='g')
                        compensate_plot(axL[0], -AlphaF,  xvals, line_from_cube,c='r')

                        verts( minmax(-AlphaF,xvals,line_from_cube))
                        verts( minmax(-AlphaF,xvals,this_c))
                        verts( minmax(-AlphaF,xvals,this_l))
                        ylim=verts.minmax
                        ylim = [1e-7,1e-3]

                        dt.axbonk(axL[0],xscale='log',yscale='log', ylim=ylim, xlim=[lcentf.min(),lcentf.max()])

                        axL[1].scatter( popt_1[1], slope)
                        dt.axbonk(axL[1],xlabel=r'$\alpha_1$',ylabel=r'$\alpha_3$', xlim=a1_ext.minmax,ylim=a3_ext.minmax)


                        axL[2].scatter( chi_l, chi_c)
                        dt.axbonk(axL[2],xscale='log',yscale='log', xlim=chi_l_ext.minmax,ylim=chi_c_ext.minmax, xlabel=r'$\chi^2_1$',ylabel=r'$\chi^3_2$')

                        axL[3].scatter( c3/b3, d3/b3)
                        dt.axbonk(axL[3], xlabel=r'$c/b$', ylabel=r'$d/b$', xlim = c_to_b.minmax, ylim=d_to_b.minmax)

                        axL[4].scatter( a3,b3)
                        dt.axbonk(axL[4],xlabel='a',ylabel='b',xlim=a_ext.minmax,ylim=b_ext.minmax)

                        counter2+=1
                        fig.savefig('%s/quad_%s_%s_%s_%02d.png'%(plotdir,field,sim,LOS, counter2))
                        


                plt.close('all')
                import three_way_bean
                plt.clf()
                fig4, ax_alpha, ax_top, ax_right = three_way_bean.three_way_bean()

                bins = np.mgrid[-5:0:16j]
                hist, bins, visthings=ax_right.hist( slope_list_3, histtype='step', orientation='horizontal',bins=bins)

                norm = mpl.colors.Normalize(vmin=min(color_3), vmax=max(color_3))
                ax_alpha.scatter( slope_list, slope_list_3, c=color_3,norm=norm)
                ax_alpha.plot( [-5,0],[-5,0],c=[0.5]*3)
                ax_top.hist( slope_list, histtype='step',bins=bins)
                dt.axbonk(ax_alpha, xlabel=r'$\alpha_1$',ylabel=r'$\alpha_3$', xlim=[-5,0],ylim=[-5,0])

                ax_right.set_ylim( ax_alpha.get_ylim())
                ax_top.set_xlim( ax_alpha.get_xlim())

                a_list=nar(a_list)
                b_list=nar(b_list)
                c_list=nar(c_list)
                d_list=nar(d_list)

                most_probable_bin = np.argmax(hist)
                bL = bins[most_probable_bin]
                bR = bins[most_probable_bin+1]
                SL = nar(slope_list_3)
                most_probable =  (SL >= bL)*(SL <= bR)
                a_mean=a_list.mean()
                b_mean=a_list.mean()
                most_probable *= ( np.abs( a_list - a_mean) < 3*a_list.std())

                AlphaC = SL[most_probable].mean()

                ax_alpha.axvline(AlphaC, c=[0.5]*3)
                ax_alpha.axhline(AlphaC, c=[0.5]*3)
                AL = nar(amp_list_3)
                Amp = AL[ most_probable ].mean()
                KL = nar(k_list_3)
                Kfit = KL[ most_probable ].mean()

                fig4.savefig('%s/alpha3_%s_%s_%s.png'%(plotdir,field,sim,LOS))
                plt.close(fig4)

                #
                # Slope Plot
                #

                fig2,ax2=plt.subplots(1,1, figsize=(12,8))
                x0=fit_range[0]; x1=fit_range[1]
                ok = (lcentf>=x0)*(lcentf<=x1)
                xvals = lcentf[ok]
                compensate_plot(ax2,-AlphaC, lcentf, specf, c=[0.5]*4)




                #
                # Fixed
                #

                this_l=10**line(np.log10(xvals), popt_1[0],AlphaF)
                kfit = lcentf[16]
                index = np.argmin( np.abs( xvals-kfit))
                this_l = this_l*specf[16]/this_l[index]
                chi2 = np.sum( (this_l-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, this_l, c='g',label=r'$%0.2e \pm %0.2e$'%(AlphaF,chi2))
                


                #
                # Most Probable
                #

                line_from_cube = 10**( AlphaC*(np.log10(xvals)-Kfit) + Amp)
                chi2 = np.sum( (line_from_cube-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, line_from_cube, c='r',label=r'$%0.2e \pm %0.2e$'%(AlphaC,chi2))

                compensate_plot(ax_all, 3, lcentf, specf, c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim])
                compensate_plot(ax_all, 3, xvals, this_l,  c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim], label=r'$\alpha_3=%0.2f$'%AlphaC)
                #ax_all.plot( lcentf, specf, c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim])
                #ax_all.plot( xvals, this_l,  c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim], label=r'$\alpha_3=%0.2f$'%AlphaC)
                ax_all.legend(loc=0)

                #
                # cube
                #

                #
                # long one.  Using the average indices makes for bad plots.
                #Most Probable
                longest_index=np.argmax(dof_list)
                a3 = a_list[longest_index]
                b3 = b_list[longest_index]
                c3 = c_list[longest_index]
                d3 = d_list[longest_index]
                k_flat = -c3/(3*d3)
                slope = b3 - c3**2/(3*d3)
                value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                line_from_cube = 10**( slope*(np.log10(xvals)-k_flat) + value)
                #line_from_cube = 10**( slope*(np.log10(xvals)-k_flat) + value)
                chi2 = np.sum( (line_from_cube-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, line_from_cube, c='k',label=r'$%0.2e \pm %0.2e$'%(slope,chi2))

                #ax2.scatter(10**(k_flat ), 10**(-AlphaC)*(value))
                #ax2.axvline(10**k_flat)
                #ax2.axhline(10**(-AlphaC)*10**value)


                #index = np.argmin( np.abs( xvals-Kfit))
                #cube_from_cube = 10**( AlphaC*(np.log10(xvals)-Kfit) + Amp)
                #cube_from_mean = 10**cube( np.log10(xvals), popt_3[0], popt_3[1], popt_3[2], popt_3[3])
                ####cube_from_mean = 10**cube( np.log10(xvals), mean_a, mean_b, mean_c, mean_d)
                cube_from_mean = 10**cube( np.log10(xvals), a3,b3,c3,d3)
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
                dt.axbonk( ax2, xscale='log',yscale='log', ylim=[1e-7,1e-3])#, ylim=ylim)
                fig2.savefig('%s/Net_%s_%s_%s.png'%(plotdir,field,sim,LOS))

                #ax_all[0].scatter( mean_a, mean_c**2/(3*mean_d), c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( AlphaC, -mean_c**2/(3*mean_d), c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( c_list/b_list,d_list/b_list, c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( c_list[most_probable]/b_list[most_probable],d_list[most_probable]/b_list[most_probable], c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( c_list[most_probable],d_list[most_probable], c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                ##ax_all[1].scatter( mean_c/mean_b,mean_d/mean_b, c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #bbb.append(mean_b); ccc.append(mean_c**2/(3*mean_d))
                SL2 = nar(c_list)**2/(3*nar(d_list))
                #ax_all[0].scatter( b_list[most_probable], SL2[most_probable],c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[0].scatter( b_list, SL2,c=sim_colors.color[sim],marker=sim_colors.marker[sim])
                #ax_all[1].scatter( quan3[sim]['msavg'], mean_a, c=sim_colors.color[sim],marker=sim_colors.marker[sim])


                MS.append( quan3[sim]['msavg'])
                MA.append( quan3[sim]['maavg'])

                #fig7,ax7=plt.subplots(1,1)
                #ax7.hist(a_list, histtype='step',color='r',linestyle='--')
                #ax7.hist(b_list, histtype='step',color='g',linestyle='--')
                #ax7.hist(c_list, histtype='step',color='b',linestyle='--')
                #ax7.hist(d_list, histtype='step',color='m',linestyle='--')
                #ax7.hist(a_list[most_probable],histtype='step',color='r',linestyle='-')
                #ax7.hist(b_list[most_probable],histtype='step',color='g',linestyle='-')
                #ax7.hist(c_list[most_probable],histtype='step',color='b',linestyle='-')
                #ax7.hist(d_list[most_probable],histtype='step',color='m',linestyle='-')
                #ax7.scatter( mean_a, 8, c='r')
                #ax7.scatter( mean_b, 8, c='g')
                #ax7.scatter( mean_c, 8, c='b')
                #ax7.scatter( mean_d, 8, c='m')
                #fig7.savefig('%s/coeff_%s_%s_%s.png'%(plotdir,field,sim,LOS))

        dt.axbonk(ax_all,xscale='log',yscale='log')
        fig_all.savefig('%s/Linearity_%s.png'%(plotdir,field))
if 0:
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

        #minmin=min([min(bbb),min(ccc)])
        #maxmax=max([max(bbb),max(ccc)])
        #bbb=nar(bbb)
        #pfit = np.polyfit( bbb,ccc,1)

        #ax_all[0].plot( bbb, bbb*pfit[0]+pfit[1])
        #print(bbb)
        ##ax_all[0].plot([minmin,maxmax],[minmin,maxmax],c=[0.5]*3)
        #dt.axbonk(ax_all[0],xlabel='b',ylabel=r'$c^2/(3 d)$')
        #fig_all.savefig('%s/Linearity_%s.png'%(plotdir,field))
        #print("FARTS %s"%field)

                #ax.plot( lcentf, specf, c=[0.5]*4)
                #ax.plot( xvals, 10**quad(np.log10(xvals), popt_2[0],popt_2[1],popt_2[2]), c='r')
                #ax.plot( xvals, 10**line(np.log10(xvals), popt_1[0],popt_1[1]),c='g')
                #dt.axbonk(ax,xscale='log',yscale='log')
                #fig.savefig('%s/quad_%s_%s_%s.png'%(plotdir,field,sim,LOS))

