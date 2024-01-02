
from GL import *
from queb3 import powerlaw_fit as plfit

import simulation

import bucket

def make_pdf(simlist):
    
    fig,ax=plt.subplots(1,1)
    for nsim,sim in enumerate(simlist):
        this_sim=simulation.corral[sim]



        for frame in this_sim.ann_frames:
            this_name = "%s/DD%04d.products/pdf_density.h5"%(this_sim.product_location,frame)
            if os.path.exists(this_name):
                print("exists, continue:",this_name)
                continue
            ds = this_sim.load_ds(frame)
            ad = ds.all_data()
            density = ad['density'].v.flatten()
            bins = np.geomspace(1e-3, 5e2,256)
            cbins=0.5*(bins[:-1]+bins[1:])
            print('pdf')
            pdf,bins = np.histogram(density,bins=bins, density=True)
            h5ptr = h5py.File(this_name,'w')
            h5ptr['hist']=pdf
            h5ptr['cbins']=cbins
            h5ptr.close()
            print("wrote",this_name)



