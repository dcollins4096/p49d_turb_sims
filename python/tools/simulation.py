from GL import *

if 'corral' not in dir():
    corral={}

class sim():
    def __init__(self,name=None,data_location=None,product_location=None, ms=None,ma=None,color='k',linestyle=':',marker="*",framelist=None, tdyn=None):
        self.name=name
        self.data_location=data_location
        self.product_location=product_location
        self.ms=ms
        self.ma=ma
        self.color=color
        self.linestyle=linestyle
        self.marker=marker
        self.framelist=framelist
        if tdyn is not None:
            self.tdyn=tdyn
        else:
            self.tdyn = 0.5/self.ms
        corral[self.name]=self

        self.quan3 = None
    def read_all_spectra(self):
        self.all_spectra={}
        for frame in self.framelist:
            self.all_spectra[frame]={}

            k3d, density = dt.dpy('%s/DD%04d.products/power_density.h5'%(self.product_location,frame), ['k','power'])
            k3d, Htotal  = dt.dpy('%s/DD%04d.products/power_Htotal.h5'%(self.product_location,frame), ['k','power'])
            k3d, velocity = dt.dpy('%s/DD%04d.products/power_velocity.h5'%(self.product_location,frame), ['k','power'])
            self.all_spectra[frame]['k3d']=k3d
            self.all_spectra[frame]['density']=density
            self.all_spectra[frame]['velocity']=velocity
            self.all_spectra[frame]['Htotal']=Htotal
            self.all_spectra[frame]['x']={}
            self.all_spectra[frame]['y']={}
            self.all_spectra[frame]['z']={}
            for axis in 'xyz':
                self.all_sepctra[frame][axis]={}
                k2d, ClTT = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTT'])
                k2d, ClEE = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClEE'])
                k2d, ClBB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClBB'])
                k2d, ClTE = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTE'])
                k2d, ClTB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTB'])
                k2d, ClEB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClEB'])
                if len(k2d) == len(ClTT)+1:
                    k2d = 0.5*(k2d[1:]+k2d[:-1])
                self.all_spectra[frame][axis]['ClTT']=ClTT
                self.all_spectra[frame][axis]['ClEE']=ClEE
                self.all_spectra[frame][axis]['ClBB']=ClBB
                self.all_spectra[frame][axis]['ClTE']=ClTE
                self.all_spectra[frame][axis]['ClTB']=ClTB
                self.all_spectra[frame][axis]['ClEB']=ClEB


        
    def read_avg_quan(self):
        if self.quan3 != None and False:
            print("Not reading twice", self.name)
            return

        self.quan3={}
        self.quan_time={}
        vrms=[]
        msavg=[]
        maavg=[]
        brms =[]
        for frame in self.framelist:
            fname = '%s/DD%04d.products/data%04d.AverageQuantities.h5'%(self.product_location,frame,frame)
            if not os.path.exists(fname):
                print("missing",fname)
                continue
            h5ptr=h5py.File(fname,'r')
            try:
                for field in h5ptr:
                    if field in self.quan_time:
                        self.quan_time[field] = np.concatenate([self.quan_time[field],h5ptr[field][()]])
                    else:
                        self.quan_time[field] = h5ptr[field][()]
                v2 = np.sqrt(h5ptr['vx_std'][:]**2+h5ptr['vy_std'][:]**2+h5ptr['vz_std'][:]**2)
                b_mean = np.sqrt(h5ptr['bx_avg'][:]**2+h5ptr['by_avg'][:]**2+h5ptr['bz_avg'][:]**2)
                b2  = np.sqrt(h5ptr['bx_std'][:]**2+h5ptr['by_std'][:]**2+h5ptr['bz_std'][:]**2)
                #NO 4 pi, this came straight off disk.
                ma=v2/b_mean
                vrms=np.append(vrms,v2)
                maavg=np.append(maavg,ma)
                brms=np.append(brms,b2)
            except:
                raise
            finally:
                h5ptr.close()
        self.quan_time['vrms']=vrms
        self.quan_time['brms']=brms
        self.quan_time['ma']=maavg
        self.quan3['maavg'] =np.mean(maavg)
        self.quan3['msavg'] =np.mean(vrms)

