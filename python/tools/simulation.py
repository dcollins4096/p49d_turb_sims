from GL import *
from queb3 import powerlaw_fit as plfit

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
        self.ann_frames=framelist
        if tdyn is not None:
            self.tdyn=tdyn
        else:
            self.tdyn = 0.5/self.ms
        corral[self.name]=self

        self.all_frames = self.get_all_frames()

        self.quan3 = None
        self.all_spectra=None
    def get_all_frames(self):
        all_dirs=glob.glob("%s/DD????"%(self.data_location))
        dir_nums = sorted([int( os.path.basename(s)[2:]) for s in all_dirs])
        if len(dir_nums)>0:
            if dir_nums[0]==0:
                dir_nums.pop(0)
        return dir_nums

    def fit_all_spectra(self):
        if self.all_spectra is None:
            print("run read_all_spectra")
            return
        self.slopes={}
        self.amps={}
        for frame in self.all_frames:
            self.slopes[frame]={}
            self.amps[frame]={}
            #3d fields first.
            k3d = self.all_spectra[frame]['k3d']
            for field in ['density','velocity','Htotal']:
                spec = self.all_spectra[frame][field]
                fitrange = [k3d[4],k3d[25]]
                slope,amp,res=plfit(k3d,spec,fitrange)
                self.slopes[frame][field]=slope
                self.amps[frame][field]=amp
            for axis in 'xyz':
                self.slopes[frame][axis]={}
                self.amps[frame][axis]={}
                k2d = self.all_spectra[frame][axis]['k2d']
                for field in ['ClTT','ClEE','ClBB']:
                    fitrange = [k2d[4], k2d[25]]
                    spec = self.all_spectra[frame][axis][field]
                    if spec is None:
                        continue
                    slope,amp,res=plfit(k2d,spec,fitrange)
                    self.slopes[frame][axis][field]=slope
                    self.amps[frame][axis][field]=amp

    def read_all_spectra(self):
        if self.all_spectra is not None:
            print("Spectra exists, not reading", self.name)
            return

        self.all_spectra={}
        for frame in self.all_frames:
            self.all_spectra[frame]={}

            k3d, density = dt.dpy('%s/DD%04d.products/power_density.h5'%(self.product_location,frame), ['k','power'])
            k3d, Htotal  = dt.dpy('%s/DD%04d.products/power_Htotal.h5'%(self.product_location,frame), ['k','power'])
            k3d, velocity = dt.dpy('%s/DD%04d.products/power_velocity.h5'%(self.product_location,frame), ['k','power'])
            self.all_spectra[frame]['k3d']=k3d.real
            self.all_spectra[frame]['density']=density.real
            self.all_spectra[frame]['velocity']=velocity.real
            self.all_spectra[frame]['Htotal']=Htotal.real
            self.all_spectra[frame]['x']={}
            self.all_spectra[frame]['y']={}
            self.all_spectra[frame]['z']={}
            for axis in 'xyz':
                self.all_spectra[frame][axis]={}
                k2d, ClTT = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTT'])
                k2d, ClEE = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClEE'])
                k2d, ClBB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClBB'])
                k2d, ClTE = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTE'])
                k2d, ClTB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTB'])
                k2d, ClEB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClEB'])
                if len(k2d) == len(ClTT)+1:
                    k2d = 0.5*(k2d[1:]+k2d[:-1])
                self.all_spectra[frame][axis]['k2d']=k2d
                self.all_spectra[frame][axis]['ClTT']=ClTT.real
                self.all_spectra[frame][axis]['ClEE']=ClEE.real
                self.all_spectra[frame][axis]['ClBB']=ClBB.real
                self.all_spectra[frame][axis]['ClTE']=ClTE.real
                self.all_spectra[frame][axis]['ClTB']=ClTB.real
                self.all_spectra[frame][axis]['ClEB']=ClEB.real


        
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
        for frame in self.all_frames:
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

