from GL import *
import yt
from volavg import *

def downsample_and_write(pf,outname, refine_by=4, write_hdf5=False,write_fits=True):
    print("RUN ", pf)

    cg = pf.covering_grid(0,[0.0]*3,pf.domain_dimensions)
    extra_dims = {'BxF':nar([1,0,0]),'ByF':nar([0,1,0]),'BzF':nar([0,0,1])}
    print("get fine")
    fine_density = cg['Density']
    print("    down density")
    coarse_density = volavg(fine_density, rank=3, refine_by = refine_by)
    if write_hdf5:
        fptr = h5py.File(outname,"w")
    for in_field_name in ['Density',
                          'x-velocity','y-velocity','z-velocity',
                          'Bx','By','Bz']: #['Density','x-velocity','y-velocity','z-velocity','BxF','ByF','BzF']:
        print("   down ",in_field_name)
        infield = cg[in_field_name]
        dims = infield.shape
        #out_field_name = {'MagneticField_F_1':'BxF','MagneticField_F_2':'ByF','MagneticField_F_3':'BzF',
        #                'Density':'Density','x-velocity':'x-velocity','y-velocity':'y-velocity','z-velocity':'z-velocity'}[in_field_name]
        out_field_name = in_field_name
        if out_field_name == 'BxF':
            field_to_store = na.ones( [(dims[0])/refine_by+1, dims[1]/refine_by, dims[2]/refine_by])
            for i in range(dims[0]/refine_by):
                field_to_store[i,:,:] = volavg(infield[refine_by*i,:,:], rank=2, refine_by=refine_by)
                field_to_store[dims[0]/refine_by,:,:] = field_to_store[0,:,:]
        elif out_field_name == 'ByF':
            field_to_store = na.ones( [(dims[0])/refine_by, (dims[1])/refine_by+1, dims[2]/refine_by])
            for j in range(dims[1]/refine_by):
                field_to_store[:,j,:] = volavg(infield[:,refine_by*j,:], rank=2, refine_by=refine_by)
                field_to_store[:,dims[1]/refine_by,:] = field_to_store[:,0,:]
        elif out_field_name == 'BzF':
            field_to_store = na.ones( [(dims[0])/refine_by, dims[1]/refine_by, (dims[2])/refine_by+1])
            for k in range(dims[2]/refine_by):
                field_to_store[:,:,k] = volavg(infield[:,:,refine_by*k], rank=2, refine_by=refine_by)
                field_to_store[:,:,dims[2]/refine_by] = field_to_store[:,:,0]
        elif out_field_name == 'Density':
            field_to_store = coarse_density
        elif out_field_name in ['x-velocity','y-velocity','z-velocity']:
            field_to_store = volavg(fine_density*infield,rank=3,refine_by=refine_by)/coarse_density
        else:
            field_to_store = volavg(infield,rank=3,refine_by=refine_by)
        if write_hdf5:
            fptr.create_dataset(in_field_name, data=field_to_store)
        if write_fits:
            outfile = "%s_%s.fits"%(outname,out_field_name)
            hdu = pyfits.PrimaryHDU(field_to_store)
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto(outfile,overwrite=True)
    if write_hdf5:
        fptr.close()

