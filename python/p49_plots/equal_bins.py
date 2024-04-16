
from GL import *

import simulation
reload(simulation)

def bin(sim_list):


    for sim in sim_list:
        this_sim = simulation.corral[sim]
        this_sim.load()

        #ann_frames is the list of analysis frames for this simulation.
        frames = this_sim.ann_frames[-1:]

        b = this_sim.get_field('b', frames[0], ax='y')

        fig,ax=plt.subplots(2,2)
        ax0=ax[0][0];ax1=ax[0][1]; ax2=ax[1][1]
        Nbins=10
        ax0.hist(b.flatten(),density=True,bins=Nbins)
        print(bin_edges)
        y = b.flatten()+0
        y.sort()
        nelements = y.size/Nbins
        a = np.arange(0,y.size+nelements,nelements,dtype='int')
        bin_edges = y[[tuple(a)]][0]
        db = (bin_edges[1:]-bin_edges[:-1])
        bc = 0.5*(bin_edges[1:]+bin_edges[:-1])
        print(db)
        ax2.bar(bc,1/(y.size*db),width=db)
        hist,bins,things = ax1.hist(y,bins=bin_edges, histtype='step',density=True)
        fig.savefig('%s/pdfs'%plotdir)


def parabin(sim_list):
    #
    # Parallel tool for fixed-frequency PDF.
    #


    for sim in sim_list:
        this_sim = simulation.corral[sim]
        this_sim.load()

        #ann_frames is the list of analysis frames for this simulation.
        frames = this_sim.ann_frames[-1:]

        b = this_sim.get_field('b', frames[0], ax='y')
        t = this_sim.get_field('d', frames[0], ax='y')
        e = this_sim.get_field('e', frames[0], ax='y')


        Nb = 1024
        NB = 16
        D = b.flatten()+0
        D.sort()

        #done in parallel
        #collect, broadcast Dmin and Dmax
        Dmin = D.min()
        Dmax = D.max()

        #fw histogram, on cores
        fw_bins = np.linspace(Dmin,Dmax,Nb)
        fw_bc = 0.5*(fw_bins[1:]+fw_bins[:-1])
        fw_db = (fw_bins[1:]-fw_bins[:-1])
        fw_hist, bins = np.histogram(D, bins=fw_bins)

        #
        # Broadcast node histogram to root histogram.
        # Sum up.
        #

        fw_pdf = fw_hist/(D.size*fw_db)
        H = np.cumsum(fw_hist)
        print("Norm fw",(fw_pdf*fw_db).sum())
        
        #The elements for the ends of the FF PDF
        J = np.linspace(0,D.size-1,NB,dtype='int')

        #ff histogram, serial, to check our result
        ff_bins = D[J]
        ff_bc = 0.5*(ff_bins[1:]+ff_bins[:-1])
        ff_db =     (ff_bins[1:]-ff_bins[:-1])
        ff_pdf = 1./(ff_db*(NB-1))
        print("Norm ff",(ff_pdf*ff_db).sum())
        print(D.size)

        #For each bin edge j, get the FW histogram bin that j is in.
        #Collect that bin to the kth processor.  Sort.  Get the right value.
        fw_bin_for_j=[]
        fw_offset = []
        for nj in range( len(J)):
            #get the first element for which Delta is strictly positive.
            #this logic can be on root, then the list of this_fw_bin[nj] and offset[nj] need to be sent.  
            j = J[nj]
            Delta = H-j
            this_fw_bin = np.where(Delta > 0 )[0][0]
            offset = j - (H[this_fw_bin]-fw_hist[this_fw_bin])
            fw_bin_for_j.append(this_fw_bin)
            fw_offset.append(offset)

        #
        # BROADCAST fw_bin_for_j and fw_offset to nodes
        #
        FF_edges_parallel = []
        for nj in range( len(J)):
            this_fw_bin = fw_bin_for_j[nj]
            offset = fw_offset[nj]

            
            #this is where each processor collects the contents of the k=this_fw_bin, 
            #send to the kth processor.
            bin_j_bool = (D>=fw_bins[this_fw_bin])*(D<=fw_bins[this_fw_bin+1])
            bin_j = D[bin_j_bool]

            #
            # Send bin_j to the kth_processor
            #

            #
            # On the kth processor,
            #

            bin_j.sort()
            element = bin_j[offset]
            FF_edges_parallel.append(element)
            got = np.where( D==element)[0][0]
            #print( ff_bins[nj], element)

        FF_edges_parallel = nar(FF_edges_parallel)
        db = FF_edges_parallel[1:]-FF-edges_parallel[:-]
        pdf = 1/db


        fig,axes=plt.subplots(1,1)
        ax0=axes#[0]#;ax1=axes[1]
        ax0.bar(fw_bc, fw_pdf, width=fw_db)
        ax0.bar(ff_bc, ff_pdf, width=ff_db, fill=False)
        fig.savefig('%s/parallel_pdf'%plotdir)






