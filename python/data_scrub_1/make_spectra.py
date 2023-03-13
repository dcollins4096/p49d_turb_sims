


from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(queb3)
reload(sim_colors)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)

clobber=False
simlist=sim_colors.simlist
#simlist=['5_half','5_1','5_2', '5_3']
shortprefix="time"
def make_spectra_files():
    for axes in ['x','y','z']:
        for i,sim in enumerate(simlist):
            print('SPECTRA',sim)
            frames=sim_colors.framelist[i]
            simdes=sim
            spectra_fname = "avg_spectra_%s_%s.h5"%(simdes,axes)
            frbname=""
            #sim_dir = "/archive2/kas14d/512reruns/frbs/%s"%simdes
            sim_dir = "/data/cb1/Projects/P49_EE_BB/%s"%sim
            product_dir = "/data/cb1/Projects/P49_EE_BB/Products/%s"%sim
            plot_dir = "./plots"
            gen_dir = "./plots"

            avg_clee=0
            avg_clbb=0
            avg_cleb=0
            avg_clte=0
            avg_cltb=0
            avg_cltt=0
            avg_v=0
            avg_h=0
            avg_d=0
            avg_rte=0
            avg_rtb=0
            avg_reb=0

            if 1:
                var_clee=0
                var_clbb=0
                var_cleb=0
                var_clte=0
                var_cltb=0
                var_cltt=0
                var_v=0
                var_h=0
                var_d=0
                var_s_clee=0
                var_s_clbb=0
                var_s_cleb=0
                var_s_clte=0
                var_s_cltb=0
                var_s_cltt=0
                var_s_v=0
                var_s_h=0
                var_s_d=0
            projections=[]
            longprefix='%s_%s_%s'%(simdes,shortprefix,axes)
            pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=longprefix, product_directory=product_dir)
            nplots=0
            proj=pack.read_queb(frame=frames[0],ax=axes,bin_style='dx1')
            fitrange=proj.determine_fit_range()  #something better
            if not os.path.exists(spectra_fname) or clobber:

                for frame in frames: 
                    
                #print(frame,simdes)
                    proj=pack.read_queb(frame=frame,ax=axes,bin_style='dx1')
                    if proj is None:
                        continue
                    nplots+=1

                    proj.read_spectra(frame, directory=product_dir)
                    projections.append(proj)
                    projections[-1].compute_harmonic_products()
                    #    ax.plot(proj.lcent,proj['ClEE'])
                    if verbose:
                        print('proj[ClEE]')
                        print(proj['ClEE'])
                    avg_cltt += proj['ClTT']
                    avg_clee += proj['ClEE']
                    avg_clbb += proj['ClBB']
                    avg_clte += proj['ClTE']
                    avg_cleb += proj['ClEB']
                    avg_cltb += proj['ClTB']
                    avg_v1 = proj['vspec'][1]
                    avg_v += np.abs(avg_v1[1:])
                    avg_h1 = proj['hspec'][1]
                    avg_h += np.abs(avg_h1[1:])
                    avg_d1 = proj['dspec'][1]
                    avg_d += np.abs(avg_d1[1:])
                    if 1:
                        var_cltt += proj['ClTT']**2 
                        var_clee += proj['ClEE']**2
                        var_clbb += proj['ClBB']**2
                        var_clte += proj['ClTE']**2
                        var_cleb += proj['ClEB']**2
                        var_cltb += proj['ClTB']**2
                        var_s_cltt += proj['var_s_ClTT']
                        var_s_clee += proj['var_s_ClEE']
                        var_s_clbb += proj['var_s_ClBB']
                        if np.isnan( var_s_cltt).any():
                            pdb.set_trace()
                        var_d += avg_d**2
                        var_v += avg_v**2
                        var_h += avg_h**2
                    avg_rte += proj['ClTE']/((proj['ClTT']*proj['ClEE'])**0.5)
                    avg_rtb += proj['ClTB']/((proj['ClTT']*proj['ClBB'])**0.5)
                    avg_reb += proj['ClEB']/((proj['ClEE']*proj['ClBB'])**0.5)
                    if np.isnan(avg_rtb).any():
                        pdb.set_trace()
                    fitrange=proj.determine_fit_range()  #something better
                    #        pack.make_spectra(frame=frame)
                    #       proj.read_spectra(frame)
                    #        pack.plot_many_spectra(proj,simdes,fname='%s/%s_spectra_%04d_%s.png'%(plot_dir,simdes,frame,axes))
                    #do the fit.  Not much fancy here
                slopes=proj.fit_eb_slopes(fitrange=fitrange)
                avg_clee /= nplots
                avg_clbb /= nplots
                avg_cleb /= nplots
                avg_clte /= nplots
                avg_cltb /= nplots
                avg_cltt /= nplots
                avg_v    /= nplots
                avg_h    /= nplots
                avg_d    /= nplots
                avg_rte  /= nplots
                avg_rtb  /= nplots
                avg_reb  /= nplots
                if 1:
                    std_cltt = np.sqrt(var_cltt/nplots- avg_cltt**2)
                    std_clee = np.sqrt(var_clee/nplots- avg_clee**2)
                    std_clbb = np.sqrt(var_clbb/nplots- avg_clbb**2)
                    std_clte = np.sqrt(var_clte/nplots- avg_clte**2)
                    std_cleb = np.sqrt(var_cleb/nplots- avg_cleb**2)
                    std_cltb = np.sqrt(var_cltb/nplots- avg_cltb**2)
                    std_d = np.sqrt(var_d/nplots- avg_d**2)
                    std_v = np.sqrt(var_v/nplots- avg_v**2)
                    std_h = np.sqrt(var_h/nplots- avg_h**2)
                    std_s_cltt = np.sqrt(var_s_cltt/nplots- avg_cltt**2)
                    std_s_clee = np.sqrt(var_s_clee/nplots- avg_clee**2)
                    std_s_clbb = np.sqrt(var_s_clbb/nplots- avg_clbb**2)


                if np.isnan(avg_rtb).any():
                    pdb.set_trace()
                fptr = h5py.File(spectra_fname,"w")
                try:
                    fptr.create_dataset("avg_cltt",data = avg_cltt)
                    fptr.create_dataset("avg_clee",data = avg_clee)
                    fptr.create_dataset("avg_clbb",data = avg_clbb)
                    fptr.create_dataset("avg_clte",data = avg_clte)
                    fptr.create_dataset("avg_cleb",data = avg_cleb)
                    fptr.create_dataset("avg_cltb",data = avg_cltb)
                    fptr.create_dataset("avg_v"   ,data = avg_v   )
                    fptr.create_dataset("avg_h"   ,data = avg_h   )
                    fptr.create_dataset("avg_d"   ,data = avg_d   )
                    fptr.create_dataset("avg_rte" ,data = avg_rte )
                    fptr.create_dataset("avg_rtb" ,data = avg_rtb )
                    fptr.create_dataset("avg_reb" ,data = avg_reb )
                    if 1:
                        fptr.create_dataset("std_cltt",data = std_cltt)
                        fptr.create_dataset("std_clee",data = std_clee)
                        fptr.create_dataset("std_clbb",data = std_clbb)
                        fptr.create_dataset("std_clte",data = std_clte)
                        fptr.create_dataset("std_cleb",data = std_cleb)
                        fptr.create_dataset("std_cltb",data = std_cltb)
                        fptr.create_dataset("std_v"   ,data = std_v   )
                        fptr.create_dataset("std_h"   ,data = std_h   )
                        fptr.create_dataset("std_d"   ,data = std_d   )
                        fptr.create_dataset("std_s_cltt",data = std_s_cltt)
                        fptr.create_dataset("std_s_clee",data = std_s_clee)
                        fptr.create_dataset("std_s_clbb",data = std_s_clbb)

                finally:
                    fptr.close()

make_spectra_files()
