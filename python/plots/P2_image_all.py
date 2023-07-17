from GL import *
import simulation

def image(sim_name):
    this_sim=simulation.corral[sim_name]

    for frame in this_sim.all_frames:
        print("Image %s %d"%(sim_name,frame))
        fig,axes=plt.subplots(3,4,figsize=(12,8))
        directory = "%s/DD%04d.products"%(this_sim.product_location,frame)
        for nax,axis in enumerate('xyz'):
            for nf,field in enumerate(['density_','magnetic_field_strength_','E','B']):
                myax=axes[nax][nf]
                frb_name = "%s/DD%04d_%s%s.fits"%(directory,frame,field,axis)
                array=pyfits.open(frb_name)[0].data
                if nf == 0:
                    array=np.log10(array)
                plot=myax.imshow(array)
                fig.colorbar(plot,ax=myax)
                myax.set(xticks=[],yticks=[])
        for nax,axis in enumerate('xyz'):
            myax=axes[nax][0].set(ylabel='Axis %s'%axis)
        for nf,field in enumerate(['density','Hmag','E','B']):
            myax=axes[0][nf]
            myax.set_title(field)
            #for nf,field in enumerate(['density_','E','B']):

        outname='%s/image_%s_n%04d.png'%(dl.plotdir,sim_name,frame)
        fig.tight_layout()
        fig.savefig(outname)




