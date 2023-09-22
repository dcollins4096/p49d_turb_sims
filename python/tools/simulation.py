from GL import *
from queb3 import powerlaw_fit as plfit

if 'corral' not in dir():
    corral={}

def set_colors(the_sim,cmap_name='jet'):
    """sets the colors on a simulation based on the mean mach number.  
    This is not a member function so it can be easily altered if necessary."""

    norm = mpl.colors.Normalize(0,7)
    cmap_function = mpl.cm.get_cmap(cmap_name)
    the_sim.color=cmap_function(norm(the_sim.Ms_mean))
    marker_size = the_sim.Ma_mean**2
    the_sim.marker='o'
    the_sim.marker_size=marker_size


class sim():
    def __init__(self,name=None,data_location=None,product_location=None, ms=None,ma=None,color='k',linestyle=':',marker="*",framelist=None, tdyn=None):
        self.name=name
        self.data_location=data_location
        self.product_location=product_location
        #self.ms=ms
        #self.ma=ma
        self.Ms_nom=ms
        self.Ma_nom=ma
        self.B_nom = ms*root4pi/ma
        self.color=color
        self.linestyle=linestyle
        self.marker=marker
        self.ann_frames=framelist
        if tdyn is not None:
            self.tdyn=tdyn
        else:
            self.tdyn = 0.5/self.Ms_nom
        corral[self.name]=self

        self.all_frames = self.get_all_frames()
        self.ann_frame_mask = nar([frame in self.ann_frames for frame in self.all_frames])

        self.quan_time = None
        self.all_spectra=None
        self.avg_spectra=None
        self.pdfs=None
        self.slopes=None
        self.slopesA=None
        self.ampsA=None
    def get_fitrange(self,xvals):
        f0 = xvals[4]
        f1 = xvals[25]
        return [f0,f1]
    def load_ds(self,frame):
        ds_name = "%s/DD%04d/data%04d"%(self.data_location,frame,frame)
        ds=yt.load(ds_name)
        return ds
    def return_queb3_package(self):
        package=queb3.simulation_package( directory=self.data_location,frames=self.all_frames,
                                                                         product_directory=self.product_location, simname=self.name)
        return package
    def get_all_frames(self):
        all_dirs=glob.glob("%s/DD????"%(self.data_location))
        dir_nums = sorted([int( os.path.basename(s)[2:]) for s in all_dirs])
        if len(dir_nums)>0:
            if dir_nums[0]==0:
                dir_nums.pop(0)
        return dir_nums
    def load(self):
        self.read_avg_quan()
        self.read_all_spectra()
        self.compute_pearson()
        self.average_spectra()
        self.fit_all_spectra()

    def fit_all_spectra(self):
        if self.all_spectra is None:
            print("run read_all_spectra")
            return
        if self.slopes is not None:
            return
        self.slopes={}
        self.amps={}
        self.slopesL = defaultdict(list)
        self.ampsL   = defaultdict(list)

        for frame in self.ann_frames:
            self.slopes[frame]={}
            self.amps[frame]={}
            #3d fields first.
            k3d = self.all_spectra[frame]['k3d']
            k2d = self.all_spectra[frame]['k2d']
            for field in self.products_positive:
                if field in self.products_3d:
                    xvals = k3d
                else:
                    xvals = k2d

                spec = self.all_spectra[frame][field]
                fitrange = self.get_fitrange(xvals)
                slope,amp,res=plfit(xvals,spec,fitrange)
                self.slopes[frame][field]=slope
                self.amps[frame][field]=amp
                self.slopesL[field].append(slope)
                self.ampsL[field].append(amp)

        if self.avg_spectra is not None:
            self.ampsA={}
            self.slopesA={}
            for field in self.products_positive:
                if field in self.products_3d:
                    xvals = k3d
                else:
                    xvals = k2d
                spec = self.avg_spectra[field]
                fitrange = self.get_fitrange(xvals)
                #fitrange = [xvals[4],xvals[25]]
                slope,amp,res=plfit(xvals,spec,fitrange)
                self.slopesA[field]=slope
                self.ampsA[field]=amp


    def compute_pearson(self):
        for frame in self.ann_frames:
            for LOS in 'xyz':
                for X,Y in [['T','E'],['T','B'],['E','B']]:
                    ClXY = self.all_spectra[frame]["Cl%s%s%s"%(X,Y,LOS)]
                    ClXX = self.all_spectra[frame]["Cl%s%s%s"%(X,X,LOS)]
                    ClYY = self.all_spectra[frame]["Cl%s%s%s"%(Y,Y,LOS)]
                    self.all_spectra[frame]['r_%s%s%s'%(X,Y,LOS)]=ClXY/np.sqrt(ClXX*ClYY)

    def average_spectra(self):
        if self.avg_spectra is not None:
            return
        self.avg_spectra={}
        for field in list(self.products)+['k2d','k3d']+self.pearson:
            self.avg_spectra[field]=0
            nframes=len(self.ann_frames)
            for frame in self.ann_frames:
                self.avg_spectra[field]=self.all_spectra[frame][field] + self.avg_spectra[field]
            self.avg_spectra[field]/=nframes

    def read_all_spectra(self):
        if self.all_spectra is not None:
            return

        self.products=['density','velocity','magnetic',
                       'ClTTx','ClTTy','ClTTz',
                       'ClEEx','ClEEy','ClEEz',
                       'ClBBx','ClBBy','ClBBz',
                       'ClTEx','ClTEy','ClTEz',
                       'ClTBx','ClTBy','ClTBz',
                       'ClEBx','ClEBy','ClEBz']
        self.products_3d = ['density','velocity','magnetic']
        self.products_positive = ['density','velocity','magnetic',
                                  'ClTTx','ClTTy','ClTTz',
                                  'ClEEx','ClEEy','ClEEz',
                                  'ClBBx','ClBBy','ClBBz']
        self.pearson=['r_TEx','r_TBx','r_EBx',
                      'r_TEy','r_TBy','r_EBy',
                      'r_TEz','r_TBz','r_EBz']

        self.all_spectra={}
        for frame in self.all_frames:
            self.all_spectra[frame]={}

            k3d, density = dt.dpy('%s/DD%04d.products/avg_power_density.h5'%(self.product_location,frame), ['k','power'])
            k3d, magnetic  = dt.dpy('%s/DD%04d.products/avg_power_magnetic.h5'%(self.product_location,frame), ['k','power'])
            k3d, velocity = dt.dpy('%s/DD%04d.products/avg_power_velocity.h5'%(self.product_location,frame), ['k','power'])
            #k3d, density = dt.dpy('%s/DD%04d.products/power_density.h5'%(self.product_location,frame), ['k','power'])
            #k3d, magnetic  = dt.dpy('%s/DD%04d.products/power_magnetic.h5'%(self.product_location,frame), ['k','power'])
            #k3d, velocity = dt.dpy('%s/DD%04d.products/power_velocity.h5'%(self.product_location,frame), ['k','power'])
            self.all_spectra[frame]['k3d']=k3d.real
            self.all_spectra[frame]['density']=density.real
            self.all_spectra[frame]['velocity']=velocity.real
            self.all_spectra[frame]['magnetic']=magnetic.real
            for axis in 'xyz':
                k2d, ClTT = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTT'])
                k2d, ClEE = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClEE'])
                k2d, ClBB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClBB'])
                k2d, ClTE = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTE'])
                k2d, ClTB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClTB'])
                k2d, ClEB = dt.dpy('%s/DD%04d.products/DD%04d_power2d%s.h5'%(self.product_location,frame,frame, axis), ['k','ClEB'])
                if len(k2d) == len(ClTT)+1:
                    k2d = 0.5*(k2d[1:]+k2d[:-1])
                self.all_spectra[frame]['k2d']=k2d
                self.all_spectra[frame]['ClTT'+axis]=ClTT.real
                self.all_spectra[frame]['ClEE'+axis]=ClEE.real
                self.all_spectra[frame]['ClBB'+axis]=ClBB.real
                self.all_spectra[frame]['ClTE'+axis]=ClTE.real
                self.all_spectra[frame]['ClTB'+axis]=ClTB.real
                self.all_spectra[frame]['ClEB'+axis]=ClEB.real


        
    def read_avg_quan(self):
        if self.quan_time != None:
            return

        self.quan_mean={} # the mean over the analysis frames.
        self.quan_time={}
        vrms=[]
        msavg=[]
        maavg=[]
        brms =[]
        frames=[]
        UNITS = root4pi
        UNITS = 1
        for frame in self.all_frames:
            frames.append(frame)
            fname = '%s/DD%04d.products/data%04d.AverageQuantities.h5'%(self.product_location,frame,frame)
            if not os.path.exists(fname):
                print("missing",fname)
                continue
            h5ptr=h5py.File(fname,'r')
            try:
                for field in h5ptr:
                    if field.startswith('b') or field.startswith('alf'):
                        U = UNITS
                    else:
                        U=1
                    if field in self.quan_time:
                        self.quan_time[field] = np.concatenate([self.quan_time[field],U*h5ptr[field][()]])
                    else:
                        self.quan_time[field] = U*h5ptr[field][()]
                v2 = np.sqrt(h5ptr['vx_std'][:]**2+h5ptr['vy_std'][:]**2+h5ptr['vz_std'][:]**2)
                b_mean = UNITS*np.sqrt(h5ptr['bx_avg'][:]**2+h5ptr['by_avg'][:]**2+h5ptr['bz_avg'][:]**2)
                b2  = UNITS*np.sqrt(h5ptr['bx_std'][:]**2+h5ptr['by_std'][:]**2+h5ptr['bz_std'][:]**2)
                #NO 4 pi, this came straight off disk.
                ma=v2/b_mean
                vrms=np.append(vrms,v2)
                maavg=np.append(maavg,ma)
                brms=np.append(brms,b2)
            except:
                raise
            finally:
                h5ptr.close()
        self.quan_time['frames']=frames
        self.quan_time['vrms']=vrms
        self.quan_time['brms']=brms
        self.quan_time['ma']=maavg
        for q in self.quan_time:
            self.quan_time[q]=nar(self.quan_time[q])
            if q[-3:]=='avg':
                self.quan_mean[q] = self.quan_time[q][self.ann_frame_mask].mean()
            elif q[-3:]=='std':
                self.quan_mean[q] = np.sqrt((self.quan_time[q][self.ann_frame_mask]**2).mean())

        self.quan_mean['maavg'] =np.mean(maavg)
        self.quan_mean['msavg'] =np.mean(vrms)
        if 1:
            self.Ma_mean = self.quan_time['ma'][self.ann_frame_mask].mean()
            self.Ms_mean = self.quan_time['vrms'][self.ann_frame_mask].mean()
        set_colors(self)

    def read_pdfs(self,fields, pdf_prefix='pdf_scaled'):
        if self.pdfs == None:
            print('read pdfs')

        self.pdfs = {}
        self.avg_pdf={}
        for field in fields:
            if field not in self.pdfs:
                self.pdfs[field]={}
            avg_pdf=0
            Npdf = 0
            for frame in self.all_frames:
                h5name = "%s/DD%04d.products/%s_%s.h5"%(self.product_location, frame,pdf_prefix, field)
                do_continue= not os.path.exists(h5name)
                if do_continue:
                    continue
                h5ptr = h5py.File(h5name,'r')
                try:
                    hist = h5ptr['hist'][()]
                    cbins= h5ptr['cbins'][()]
                except:
                    raise
                finally:
                    h5ptr.close()
                    self.pdfs[field][frame]={'hist':hist,'cbins':cbins}
                    avg_pdf = hist + avg_pdf
                    Npdf+=1
            self.avg_pdf[field]=avg_pdf/Npdf




