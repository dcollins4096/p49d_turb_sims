

from GL import *

import simulation

def slope_time(simlist, THRESHOLD=5):
    plt.close('all')
    for this_simname in simlist:
        this_sim = simulation.corral[this_simname]
        #this_sim.read_all_spectra()
        #this_sim.fit_all_spectra()
        this_sim.load()


        total_columns = 3
        total_rows = len(this_sim.products_positive)//total_columns
        fig,axes=plt.subplots(total_rows,total_columns, figsize=(8,12))
        times = this_sim.quan_time['time']
        firsts=np.zeros(len(this_sim.products_positive))
        for nf,field in enumerate(this_sim.products_positive):
            ncol = nf%total_columns
            nrow = nf//total_columns
            thax=axes[nrow][ncol]
            #thax.plot( times, this_sim.slopesL[field])
            slopes=nar(this_sim.slopesL[field])
            nslope=len(slopes)
            half_avg = np.mean(slopes[nslope//2:])
            half_std = np.std(slopes[nslope//2:])
            renorm = np.abs((slopes-half_avg)/half_std)
            thax.axhline(2)
            thax.axhline(5)
            first_good_frame = np.where( renorm <THRESHOLD)[0][0]
            firsts[nf]=this_sim.quan_time['frames'][first_good_frame]
            thax.plot(  slopes, marker='*')
            thax.plot(  renorm, marker='*')
            thax.set(title=field,xlabel='t')
        print("First Good Frame %s %d"%(this_sim.name, firsts.max()))
        outname = "%s/slope_time_%s"%(dl.plotdir,this_simname)
        fig.savefig(outname)
        print(outname)








