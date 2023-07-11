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
    def read_avg_quan(self):
        if self.quan3 != None and False:
            print("Not reading twice", self.name)
            return

        self.quan3={}
        self.quan_time={}
        v2avg=[]
        msavg=[]
        maavg=[]
        b2avg =[]
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
                ma=v2/b_mean
                v2avg=np.append(v2avg,v2)
                maavg=np.append(maavg,ma)
                b2avg=np.append(b2avg,b2)
            except:
                raise
            finally:
                h5ptr.close()
        self.quan_time['vrms']=v2avg
        self.quan_time['brms']=np.sqrt(b2avg)
        self.quan_time['ma']=maavg
        self.quan3['maavg'] =np.mean(maavg)
        self.quan3['msavg'] =np.mean(v2avg)

