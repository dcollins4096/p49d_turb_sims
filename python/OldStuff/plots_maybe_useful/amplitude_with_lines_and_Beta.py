
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

#
# Read primitives
# 

if 'quand' not in dir():
    quand = {}
for i, sim in enumerate(sim_colors.simlist):
    if sim in quand:
        continue
    print("=== %s ==="%sim)
    ol = "%s/%s/OutputLog"%(sim_dir,sim)
    quand[sim] = gaq.get_quantities_and_rms(ol)

#
# Read spectra
#

plotdir = "//home/dccollins/PigPen"
if 'clobber' not in dir():
    clobber=False
if 'spectra_dict' not in dir() or clobber==True:
    spectra_dict={}
    shortprefix='time'
    for axes in ['x','y','z']:
        print('axes',axes)
        spectra_dict[axes]={}
        for i, sim in enumerate(sim_colors.simlist):
            spectra_dict[axes][sim]=slopes(sim) 

            sim_dir = "/data/cb1/Projects/P49_EE_BB/512_frbs/%s"%sim
            frb_name = ""
            longprefix='%s_%s_%s'%(sim,shortprefix,axes)
            spectra_fname = "avg_spectra_%s_%s.h5"%(sim,axes)
            pack = queb3.simulation_package( directory=sim_dir,frames=sim_colors.frames[sim],prefix=longprefix,
                                           frbname=frb_name)
            proj=pack.read_queb(frame=sim_colors.frames[sim][0],ax=axes,bin_style='dx1')
            fitrange=proj.determine_fit_range()  #something better
            if os.path.exists(spectra_fname):
                try:
                    fptr=h5py.File(spectra_fname,'r')
                    for field in fptr:
                        spectra_dict[axes][sim].spectra[field]=fptr[field][()]
                except:
                    raise
                finally:
                    fptr.close()
            spectra_dict[axes][sim].fit(proj.lcent,fitrange)

#
# Read analysis set.
#

if 'quan3' not in dir():
    quan3={}
    for ns, sim in enumerate(sim_colors.simlist):
        quan3[sim]={}
        v2avg=[]
        msavg=[]
        maavg=[]
        for frame in sim_colors.framelist[ns]:
            fname = '%s/%s/DD%04d/data%04d.AverageQuantities.h5'%(sim_colors.cloudbreak_base,sim,frame,frame)
            if not os.path.exists(fname):
                print(fname)
                continue
            quan4=h5py.File(fname,'r')
            v2 = np.sqrt(quan4['vx_std'][:]**2+quan4['vy_std'][:]**2+quan4['vz_std'][:]**2)
            bt = np.sqrt(quan4['bx_avg'][:]**2+quan4['by_avg'][:]**2+quan4['bz_avg'][:]**2)
            ma=v2/bt
            v2avg=np.append(v2avg,v2)
            maavg=np.append(maavg,ma)
            print(np.mean(ma), np.mean(maavg))
        quan3[sim]['maavg'] =np.mean(maavg)
        quan3[sim]['msavg'] =np.mean(v2avg)


#
# Primitives Lines
#
AAA=defaultdict(list)
MMMa=defaultdict(list)
MMMs=defaultdict(list)
matrix_d = np.zeros([4,3]) #ms,ma
matrix_v = np.zeros([4,3]) #ms,ma
matrix_h = np.zeros([4,3]) #ms,ma
the_x = np.zeros_like(matrix)
the_y = np.zeros_like(matrix)
if 1:
    plt.close('all')
    fig,ax = plt.subplots(3,2, figsize=(8,4))#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axis = 'y'
    simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
    sublist=nar([["half_half","half_1","half_2"],["1_half","1_1","1_2"],["2_half","2_1","2_2"],["3_half","3_1","3_2"]]).transpose()

    for nsub,sub  in enumerate(sublist):

        # Nominal; measure with analysis set; measured with all frames (wonky)
        
        if 1:
            this_Ms =[ quan3[sim]['msavg'].mean() for sim in sub]
            this_Ma =[ quan3[sim]['maavg'].mean() for sim in sub]
            measured_or_nominal = 'measured3'
        elif 0:
            this_Ms =[ quand[sim]['Ms'].mean() for sim in sub]
            this_Ma =[ quand[sim]['Ma'].mean() for sim in sub]
            measured_or_nominal = 'measured'
        else:
            this_Ms = [sim_colors.Ms[sim]]*len(sub)
            this_Ma = [sim_colors.Ma[sim]]*len(sub)
            measured_or_nominal = 'nominal'


        #
        # slopes or amplitudes
        #
        do_log=False
        if 0:
            do_log=True
            avg_v  = [spectra_dict[axis][sim].amps['avg_v'] for sim in sub]
            avg_d  = [spectra_dict[axis][sim].amps['avg_d'] for sim in sub]
            avg_h  = [spectra_dict[axis][sim].amps['avg_h'] for sim in sub]
            aos = 'amps'
            Aalpha = 'A'
        else:
            avg_v  = [spectra_dict[axis][sim].slopes['avg_v'] for sim in sub]
            avg_d  = [spectra_dict[axis][sim].slopes['avg_d'] for sim in sub]
            avg_h  = [spectra_dict[axis][sim].slopes['avg_h'] for sim in sub]
            aos = 'slopes'
            Aalpha = "\\alpha"


        MeanBeta = nar([(quand[sim]['Btot']**2/(8*np.pi)).mean() for sim in sub])
        MeanB = nar([(quand[sim]['Btot']).mean() for sim in sub])
        outname = '%s/prim_mach_%s_%s.pdf'%(plotdir, aos, measured_or_nominal)
        colors =[sim_colors.color[sim] for sim in sub]
        print(colors)
        markers=[sim_colors.marker[sim] for sim in sub]

        #Its dumb as hell that marker isn't a list
        #and I have to make this loop
        for i in range(len(markers)):
            ax[0][0].scatter(this_Ms[i], avg_d[i], c=colors[i],marker=markers[i])
            ax[1][0].scatter(this_Ms[i], avg_v[i], c=colors[i],marker=markers[i])
            ax[2][0].scatter(this_Ms[i], avg_h[i], c=colors[i],marker=markers[i])

            ax[0][1].scatter(this_Ma[i], avg_d[i], c=colors[i],marker=markers[i])
            ax[1][1].scatter(this_Ma[i], avg_v[i], c=colors[i],marker=markers[i])
            ax[2][1].scatter(this_Ma[i], avg_h[i], c=colors[i],marker=markers[i])

            #ax[0][2].scatter(  MeanBeta[i], avg_d[i], c=colors[i],marker=markers[i])
            #ax[1][2].scatter(  MeanBeta[i], avg_v[i], c=colors[i],marker=markers[i])
            #ax[2][2].scatter(  MeanBeta[i], avg_h[i], c=colors[i],marker=markers[i])

        fits_v = np.polyfit(MeanBeta, avg_v, 1)
        fits_d = np.polyfit(MeanBeta, avg_d, 1)
        fits_h = np.polyfit(MeanBeta, avg_h, 1)
        #ax[0][2].plot(MeanB, fits_v[0]*MeanB+fits_v[1], c=colors[0])
        #ax[1][2].plot(MeanB, fits_d[0]*MeanB+fits_d[1], c=colors[0])
        #ax[2][2].plot(MeanB, fits_h[0]*MeanB+fits_h[1], c=colors[0])

        matrix_d[:,nsub]=avg_d
        matrix_v[:,nsub]=avg_v
        matrix_h[:,nsub]=avg_h
        the_x[:,nsub]=MeanBeta
        the_y[:,nsub]=nar(this_Ms)


        ax[0][0].set_ylabel(r'$%s_\rho$'%Aalpha)
        ax[1][0].set_ylabel(r'$%s_v$'%Aalpha)
        ax[2][0].set_ylabel(r'$%s_H$'%Aalpha)
        ax[2][0].set_xlabel(r'$M_{\rm{s}}$')
        ax[2][1].set_xlabel(r'$M_{\rm{A}}$')
        ax[2][2].set_xlabel(r'$\beta$')

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
        print(f)
        amps = f[0] + f[1]*myx + f[2]*myy+f[3]*myx*myy
        chi2 = ( (myval-amps)**2/amps).sum()
        print("chi", c,chi2)
        aax.scatter(myx, amps, c=c,s=0.2)
    do_morestuff( ax[0][2], the_x, the_y**2, matrix_v,c='m')
    do_morestuff( ax[1][2], the_x, the_y**2, matrix_d,c='m')
    do_morestuff( ax[2][2], the_x, the_y**2, matrix_h,c='m')

    for a in [ax[0][1],ax[1][1],ax[2][1]]:
        a.set_yticks([])
    for a in [ax[0][-1],ax[1][-1],ax[2][-1]]:
        a.yaxis.tick_right()
    if do_log:
        for a in ax.flatten():
            a.set_yscale('log')

    fig.savefig(outname)
    print(outname)

if 0:
    nominal_Ma =nar([sim_colors.Ma[sim] for sim in sub])
    nominal_Ms =nar([sim_colors.Ms[sim] for sim in sub])
    ax[0][3].scatter(nominal_Ma.mean(), fits_v[0], c=colors[0])
    AAA['v'].append(fits_v[0])
    MMMa['v'].append(nominal_Ma.mean())
    MMMa['v']=nar(MMMa['v'])
    AAA['v']=nar(AAA['v'])
    pfit = np.polyfit(MMMa['v'],AAA['v'],1)
    ax[0][3].plot(MMMa['v'],pfit[0]*MMMa['v']+pfit[1],c='k')

