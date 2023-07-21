from GL import *
import yt

import simulation


def make_magnetic_pdf(simlist,frames=None, renormalize=True, pdf_prefix='pdf_scaled', clobber=False):
    fields = ['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength']
    if renormalize:
        raw_or_renorm='renorm'
    else:
        raw_or_renorm='raw'

    for sim in simlist:
        this_sim = simulation.corral[sim]
        this_sim.load()

        if frames == 'last':
            framelist = this_sim.all_frames[-1:]
        elif frames == 'last_two':
            framelist = this_sim.all_frames[-2:]
        elif frames == 'all':
            framelist = this_sim.all_frames

        NBINS=256
        if renormalize:
            def puller(arr):
                return arr[this_sim.frame_mask].mean()
            def puller2(arr):
                return np.sqrt( (arr[this_sim.frame_mask]**2).mean())
            UNITS = root4pi
            bx,by,bz =    [UNITS*puller(this_sim.quan_time['b%s_avg'%s]) for s in 'xyz']
            sbx,sby,sbz = [UNITS*puller2(this_sim.quan_time['b%s_std'%s]) for s in 'xyz']
            vx,vy,vz =    [puller(this_sim.quan_time['v%s_std'%s]) for s in 'xyz']
            mu_i = {'magnetic_field_x':bx,'magnetic_field_y':by,'magnetic_field_z':bz}
            sig_i = {'magnetic_field_x':sbx,'magnetic_field_y':sby,'magnetic_field_z':sbz}
            v3d = np.sqrt(vx*vx+vy*vy+vz*vz)
            b3d = np.sqrt(sbx**2+sby**2+sbz**2)
            mu_i = {'magnetic_field_x':bx,'magnetic_field_y':by,'magnetic_field_z':bz, 'magnetic_field_strength':0}
            sig_i = {'magnetic_field_x':sbx,'magnetic_field_y':sby,'magnetic_field_z':sbz,'magnetic_field_strength':b3d}
            LOW=-10
            HIGH=10
            all_bins = {'magnetic_field_x':np.linspace(LOW,HIGH,NBINS),
                    'magnetic_field_y':np.linspace(LOW,HIGH,NBINS),
                    'magnetic_field_z':np.linspace(LOW,HIGH,NBINS),
                    'magnetic_field_strength':np.linspace(0,2*HIGH,NBINS)}

        for frame in framelist:
            ds_name = "%s/DD%04d/data%04d"%(this_sim.data_location,frame,frame)
            ds = None
            for field in fields:
                h5name = "%s/DD%04d.products/%s_%s.h5"%(this_sim.product_location,frame,pdf_prefix, field)
                if os.path.exists(h5name) and not clobber:
                    print("File exists. Skip.", h5name)
                    continue

                if ds is None:
                    ds = yt.load(ds_name)
                    ad = ds.all_data()
                print("PDF %s %s %d"%(this_sim.name,field,frame))
                F = ad[field].v
                if renormalize:
                    mu = mu_i[field]
                    sig= sig_i[field]
                    F = (F-mu)/sig
                    bins = all_bins[field]

                if not renormalize:
                    Fmin = F.min()
                    Fmax = F.max()
                    do_log=False
                    if Fmin <= 0:
                        do_log=False
                    if do_log:
                        bins = np.geomspace(Fmin,Fmax,64)
                    else:
                        bins = np.linspace(Fmin,Fmax,64)

                hist, xbins = np.histogram(F,bins=bins)

                h5ptr = h5py.File(h5name,'w')
                h5ptr['bins']=bins
                h5ptr['cbins']=0.5*(bins[1:]+bins[:-1])
                h5ptr['hist']=hist
                if renormalize:
                    h5ptr[field+"_mu"]=mu
                    h5ptr[field+"_sig"]=sig
                h5ptr.close()
                print(h5name)

