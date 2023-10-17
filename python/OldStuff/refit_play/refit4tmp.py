
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

field_list=['avg_cltt','avg_clee','avg_clbb', 'avg_d','avg_v','avg_h']
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
simlist = ['half_half','1_half','2_half','3_half']#['1_half']#['half_2', 'half_half']#['half_half']#['1_1']#['2_2'] #['half_half'] # ['half_half']#,'2_2']
simlist = ['3_half']

field_list = ['avg_clbb']
do_all_subplots=True


#read_stuff reads spectra and average quantities
#Reload to re-read, but don't do that every time.
import read_stuff as rs
#reload(read_stuff)  
#spectra_dict = rs.spectra_dict
quan3 = rs.quan3
import refit3
#reload(refit3)

spectra_dict = refit3.spectra_dict

from scipy.optimize import curve_fit

#
# Spectra plots
#

plotdir=os.environ['HOME']+"/PigPen"

from scipy.ndimage import gaussian_filter1d
def compensate_plot(ax,alpha,x,y,**kwargs):
    ax.plot( x, x**alpha*y, **kwargs)
def compensate_scatter(ax,alpha,x,y,**kwargs):
    ax.scatter( x, x**alpha*y, **kwargs)
if 1:
    #
    # Development code.  
    #
    #fill objects with slope array

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
                mah_bucket = spectra_dict[LOS][sim].slopes3[field]
                print("DO", field, sim, LOS)
                proj=rs.proj_dict[LOS][sim]
                fit_range=proj.determine_fit_range()
                x0=fit_range[0]
                x1=fit_range[1]
                lcentf = mah_bucket.lcentf
                specf =  mah_bucket.specf

                # fixed fit
                ok = (lcentf>=x0)*(lcentf<=x1)
                xvals = lcentf[ok]
                spec=specf[ok]
                popt_fixed, pcov_fixed= curve_fit( refit3.line, np.log10(xvals), np.log10(spec))
                AlphaF = popt_fixed[1]
                AlphaComp1 = AlphaF

                counter2=0
                slope_list =  nar(mah_bucket.line_slope_list)[:,1]
                slope_list_3 =  mah_bucket.cube_slope_list
                amp_list_3 =  mah_bucket.cube_amp_list

                ok_list= mah_bucket.ok_list
                color_3=[]
                a_list=mah_bucket.cube_param_list[:,0]
                b_list=mah_bucket.cube_param_list[:,1]
                c_list=mah_bucket.cube_param_list[:,2]
                d_list=mah_bucket.cube_param_list[:,3]

                dof_list=[ok.sum() for ok in mah_bucket.ok_list]
                verts=dt.extents()
                for nlist, ok in enumerate(mah_bucket.ok_list):
                    verts( refit3.minmax( -AlphaComp1, lcentf[ok], mah_bucket.line_list[nlist]))
                    verts( refit3.minmax( -AlphaComp1, lcentf[ok], mah_bucket.cube_list[nlist]))
                c_to_b=dt.extents( c_list/b_list)
                d_to_b=dt.extents( d_list/b_list )
                a3_ext=dt.extents( a_list)
                a1_ext=dt.extents( slope_list)
                chi_c_ext=dt.extents( nar([1e-3,1e3]))
                chi_l_ext=dt.extents( nar([1e-3,1e3]))
                a_ext=dt.extents( a_list)
                b_ext=dt.extents( b_list )

                #a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                #slope = b3 - c3**2/(3*d3)
                #value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                #k_flat = -c3/(3*d3)

                fig,axes = plt.subplots(2,2,figsize=(12,12))
                axL = axes.flatten()

                axL[2].plot([1e-6,200],[1e-6,200],c='k')
                axL[2].plot([1e-6,200],[0.1*1e-6,0.1*200],c='k')
                axL[2].plot([1e-6,200],[10*1e-6,10*200],c='k')

                axL[0].scatter( slope_list, slope_list_3)
                dt.axbonk(axL[0],xlabel=r'$\alpha_1$',ylabel=r'$\alpha_3$', xlim=a1_ext.minmax,ylim=a3_ext.minmax)

                axL[1].scatter( c_list/b_list, d_list/b_list)
                dt.axbonk(axL[1], xlabel=r'$c/b$', ylabel=r'$d/b$', xlim = c_to_b.minmax, ylim=d_to_b.minmax)

                axL[3].scatter( a_list, b_list)
                dt.axbonk(axL[3],xlabel='a',ylabel='b',xlim=a_ext.minmax,ylim=b_ext.minmax)

                axL[2].scatter( mah_bucket.chi_l_list, mah_bucket.chi_c_list)
                dt.axbonk(axL[2],xscale='log',yscale='log', xlim=chi_l_ext.minmax,ylim=chi_c_ext.minmax, xlabel=r'$\chi^2_1$',ylabel=r'$\chi^3_2$')



                fig.savefig('%s/quad_%s_%s_%s.png'%(plotdir,field,sim,LOS))

                fig_sub,ax_sub_arr=plt.subplots(1,2,figsize=(12,8))
                ax_sub=ax_sub_arr[0]
                ax_sub_1=ax_sub_arr[1]

                #
                # Do a bunch of plots
                #
                if do_all_subplots:
                    for nlist,ok in enumerate(mah_bucket.ok_list):
                        #if counter2 > 4:
                        #    continue

                        xvals = mah_bucket.lcentf[ok]
                        xvals = mah_bucket.lcentf
                        popt_1 = mah_bucket.line_slope_list[nlist]
                        popt_3 = mah_bucket.cube_param_list[nlist]

                        a3,b3,c3,d3=popt_3 #a + b k + c k^2 + d k ^3
                        slope = b3 - c3**2/(3*d3)
                        value = a3 - b3*c3/(3*d3) + c3**3/d3**2*(1./9-1./27)
                        k_flat = -c3/(3*d3)
                        if 10**k_flat <= x0 or 10**k_flat >= x1:
                            continue
                        print("x0 %0.4f kf %0.4f"%(x0,10**k_flat))
                        #slope2 = mah_bucket.cube_slope_list[nlist]
                        #print( slope, slope2)
                        #if np.abs( slope-slope2)>1e-7:
                        #    raise

                        this_l=10**refit3.line(np.log10(xvals), popt_1[0],popt_1[1])
                        line_from_cube = 10**( slope*(np.log10(xvals)-k_flat) + value)
                        this_c = 10**refit3.cube( np.log10(xvals), popt_3[0], popt_3[1], popt_3[2], popt_3[3])

                        ax_sub.clear()
                        #axL[0].set_title(r'$\alpha_1 = %0.2f \alpha_3 = %0.2f \chi^2_1 = %0.1e \chi^2_3 = %0.1e$'%(popt_1[1], slope, chi_l, chi_c))

                        compensate_plot(ax_sub, -AlphaComp1, lcentf, specf, c=[0.5]*4)
                        compensate_plot(ax_sub, -AlphaComp1, xvals, this_c, c='b')
                        compensate_plot(ax_sub, -AlphaComp1,  xvals, this_l,c='g')
                        compensate_plot(ax_sub, -AlphaComp1,  xvals, line_from_cube,c='r')
                        compensate_scatter(ax_sub,-AlphaComp1, 10**k_flat, 10**value)



                        xlim = [lcentf.min(),lcentf.max()]
                        xlim = [proj.lcent.min(),proj.lcent.max()]
                        dt.axbonk(ax_sub,xscale='log',yscale='log', ylim=verts.minmax, xlim=xlim)

                        lX = np.log10(xvals)
                        D1 = 10**(b3+2*c3*lX+3*d3*lX**2)
                        D1 = (b3+2*c3*lX+3*d3*lX**2)
                        print(D1.max())
                        compensate_plot(ax_sub_1,-AlphaComp1, xvals, D1)
                        #ax_sub_1.plot( xvals, D1)







                        counter2+=1
                        fig_sub.savefig('%s/sub_spectrum_%s_%s_%s_%02d.png'%(plotdir,field,sim,LOS, counter2))

                        


                plt.close('all')
                import three_way_bean
                plt.clf()
                fig4, ax_alpha, ax_top, ax_right = three_way_bean.three_way_bean()

                bins = np.mgrid[-5:0:51j]
                hist, bins, visthings=ax_right.hist( slope_list_3, histtype='step', orientation='horizontal',bins=bins)

                #norm = mpl.colors.Normalize(vmin=min(color_3), vmax=max(color_3))
                ax_alpha.scatter( slope_list, slope_list_3)#, c=color_3,norm=norm)
                ax_alpha.plot( [-5,0],[-5,0],c=[0.5]*3)
                ax_top.hist( slope_list, histtype='step',bins=bins)
                dt.axbonk(ax_alpha, xlabel=r'$\alpha_1$',ylabel=r'$\alpha_3$', xlim=[-5,0],ylim=[-5,0])

                ax_right.set_ylim( ax_alpha.get_ylim())
                ax_top.set_xlim( ax_alpha.get_xlim())

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
                Amp = nar(amp_list_3)[ most_probable ].mean()
                #k_flat = -c3/(3*d3)
                KL = -c_list/(3*d_list)
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

                this_l=10**refit3.line(np.log10(xvals), popt_1[0],AlphaF)
                kfit = lcentf[16]
                index = np.argmin( np.abs( xvals-kfit))
                this_l = this_l*specf[16]/this_l[index]
                chi2 = np.sum( (this_l-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, this_l, c='g',label=r'$%0.2e \pm %0.2e$'%(AlphaF,chi2))

                #
                # Most Probable
                #

                probable_line_from_cube = 10**( AlphaC*(np.log10(xvals)-Kfit) + Amp)
                chi2 = np.sum( (probable_line_from_cube-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals, probable_line_from_cube, c='r',label=r'$%0.2e \pm %0.2e$'%(AlphaC,chi2))


                #
                # cube
                #

                #
                # long one.  Using the average indices makes for bad plots.
                #Most Probable
                longest_index=np.argmax(dof_list)
                print(longest_index)
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
                xvals_longest = lcentf[ ok_list[longest_index]]
                cube_from_mean = 10**refit3.cube( np.log10(xvals_longest), a3,b3,c3,d3)
                #ax2.scatter( xvals[index],xvals[index]**(-AlphaC)*cube_from_mean[index])
                #ax2.scatter( 10**Kfit, cube_from_mean.min())
                chi2 = -42 #np.sum( (cube_from_mean-specf[ok])**2/specf[ok])
                compensate_plot(ax2, -AlphaC, xvals_longest, cube_from_mean, c='b',label=r'$%0.2e \pm %0.2e$'%(AlphaC,chi2))

#                if field in ['avg_cltt','avg_clee','avg_clbb']:
#                    ylim = [1e-11,1e6] #whole range
#                    #ylim = [1e-3, 1e-2]
#                else:
#                    ylim = [1e-10,0.5]
                ax2.legend(loc=0)
                dt.axbonk( ax2, xscale='log',yscale='log')#, ylim=[1e-7,1e-3])#, ylim=ylim)
                fig2.savefig('%s/Net_%s_%s_%s.png'%(plotdir,field,sim,LOS))


                #
                # Summary plot
                #

                if 1:
                    compensate_plot(ax_all, 0, lcentf, specf, c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim])
                    #compensate_plot(ax_all, 3, xvals, this_l,  c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim], label=r'$\alpha_3=%0.2f$'%AlphaC)
                    compensate_plot(ax_all, 0, xvals, probable_line_from_cube,  c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim], label=r'$\alpha_3=%0.2f$'%AlphaC)
                    compensate_plot(ax_all, 0, lcentf, mah_bucket.specf_raw,  c=[0.5]*3)
                #ax_all.plot( lcentf, specf, c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim])
                #ax_all.plot( xvals, this_l,  c=sim_colors.color[sim],linestyle=sim_colors.linestyle[sim], label=r'$\alpha_3=%0.2f$'%AlphaC)
                ax_all.legend(loc=0)

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
