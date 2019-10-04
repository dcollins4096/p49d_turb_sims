from GL import *
import yt
import p49_QU2EB
import p49_fields
import cmbtools
import davetools as dt
reload(dt)
import spectra_tools as st
reload(st)
frbname="frbs"
reload(p49_QU2EB)


class queb_snapshot():
    """Container for QUEB, T (density), H(magnetic field).
    Contains fields and their transforms.
    """
    def __init__(self,Q,U,T=None,H=None,BoxSize=1.0,axis='x',frame=-1, simulation=None):
        self.Q=Q
        self.U=U
        self.T=T
        self.H=H
        self.BoxSize=BoxSize
        self.axis=axis
        self.frame=frame
        self.simulation=simulation
    def __getitem__(self,item):
        """square bracket access method"""
        #this is a kludge
        return self.__dict__[item]
    def write(self):
        xd='DD'
        frb_dir = "%s/%s"%(self.simulation.directory,p49_QU2EB.frbname)
        product_dir = "%s/DD%04d.products"%(self.simulation.directory,self.frame)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,self.frame,self.axis)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,self.frame,self.axis)
        Clf= "%s/%s%04d_Cl%s.fits"%(frb_dir,xd,self.frame,self.axis)
        hdu = pyfits.PrimaryHDU(self.E)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Ef,overwrite=True)
        hdu = pyfits.PrimaryHDU(self.B)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Bf,overwrite=True)
        np.savetxt(Clf, list(zip(self.lbins,self.ClEE,self.ClBB, self.ClTE, self.ClEB)))

    def compute_harmonic_products(self):
        """
        N = np.array(arr.shape,dtype = np.int32)
        arr = np.ascontiguousarray(arr,dtype=np.double)
        xsize = N.max() #5 * np.pi / 180*BoxSize
        size2d = np.array([xsize,xsize])
        Delta = size2d/N
        Deltal = cmbtools.Delta2l(Delta,N)
        harm = cmbtools.map2harm(arr,Delta)
        lmax = Deltal[0]*N[0]/2
        lbins = np.arange(N[0]//2+1)*Deltal[0] #np.linspace(0,lmax,N//2)
        lcent = lbins[:-1] + np.diff(lbins)/2.
        ClBB = cmbtools.harm2cl(harm,Deltal,lbins)
        output={'Delta':Delta,'Deltal':Deltal,'harm':harm,'lbins':lbins,'ClBB':ClBB,'lmax':lmax}
        return output
        """
        if not self.Q.flags['C_CONTIGUOUS']:
            self.Q = np.ascontiguousarray(Q)
        if not self.U.flags['C_CONTIGUOUS']:
            self.U = np.ascontiguousarray(U)
        if not self.T.flags['C_CONTIGUOUS']:
            self.T = np.ascontiguousarray(T)
        self.N = np.array(self.Q.shape,dtype = np.int32)

        #We arrange things so that Delta=Delta x = 1.
        #Thus Delta L = 2pi/N
        self.xsize = 64*self.BoxSize
        self.size2d = np.array([self.xsize]*2)
        self.Delta = self.size2d/self.N
        self.Deltal = 2*np.pi/(self.N*self.Delta) #cmbtools.Delta2l(self.Delta,self.N) #2 pi/(N Delta)
        self.lmax = self.Deltal[0]*self.N[0]/2
        self.lbins = np.arange(0,self.N[0]//2) *self.Deltal[0]
        self.lcent = self.lbins[:-1] + np.diff(self.lbins)/2.

        self.Qharm = cmbtools.map2harm(self.Q,self.Delta)
        self.Uharm = cmbtools.map2harm(self.U,self.Delta)

        self.Eharm, self.Bharm = cmbtools.QU2EB(self.Qharm,self.Uharm,self.Deltal)
        self.E = cmbtools.harm2map(self.Eharm,self.Delta)
        self.B = cmbtools.harm2map(self.Bharm,self.Delta)
        if self.T is not None:
            self.Tharm = cmbtools.map2harm(self.T,self.Delta)
        if self.H is not None:
            self.Hharm = cmbtools.map2harm(self.H,self.Delta)
        self.ClEE = cmbtools.harm2cl(self.Eharm,self.Deltal,self.lbins)
        self.ClBB = cmbtools.harm2cl(self.Bharm,self.Deltal,self.lbins)
        if self.Hharm is not None:
            self.ClHH= cmbtools.harm2cl(self.Hharm,self.Deltal,self.lbins)
        self.ClTT = cmbtools.harm2cl(self.Tharm,self.Deltal,self.lbins)
        self.ClTE = cmbtools.harm2clcross_samegrid(self.Tharm,self.Eharm,self.Deltal,self.lbins)
        self.ClTB = cmbtools.harm2clcross_samegrid(self.Tharm,self.Bharm,self.Deltal,self.lbins)
        self.ClEB = cmbtools.harm2clcross_samegrid(self.Eharm,self.Bharm,self.Deltal,self.lbins)


class simulation_package():
    """container for a simulation.
    Keeps track of the data location and BoxSize.
    Produces FRBs from simulation data."""
    def __init__(self, clobber=False,directory=".",frames=[], prefix="RUN", 
                 fit_range=None, BoxSize=1, plot_format='png'):
        self.stuff={}
        self.clobber=clobber
        self.plot_format=plot_format
        self.clobber=clobber
        self.directory=directory
        self.frames=frames
        self.prefix=prefix
        self.fit_range=fit_range
        self.BoxSize=BoxSize

    def EBall(self):
        for frame in self.frames:
            ds = yt.load("%s/DD%04d/data%04d"%(self.directory,frame,frame))
            p49_fields.add_QU(ds)
            self.make_frbs(frame,ds=ds)
            for axis in 'xyz':
                ts=self.read_queb(frame,axis)
                ts.compute_harmonic_products()
                ts.write()

    def make_frbs(self,frame, axes=['x','y','z'], ds=None):
        fields=[]
        for axis in axes:
          fields.append( (axis,'Q%s'%(axis))   )
          fields.append( (axis,'U%s'%(axis))   )
          fields.append( (axis,'density') )
          fields.append( (axis,'magnetic_field_strength'))

        for axis, field in fields :
            outputdir = "%s/%s/"%(self.directory,frbname)
            if not os.access(outputdir, os.F_OK):
                os.mkdir(outputdir)
            #fix names; Q and U have the name in the field, but others don't.
            if field[0] in 'QU' and field[1] in 'xyz':
                field_name = field
            else:
                field_name = field + "_"+axis
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field_name)
            if os.access(outfile, os.F_OK) and not self.clobber:
                print("FRB exists: %s"%outfile)
            else:
                print("FRB being produced: %s"%outfile)
                res = ds.parameters['TopGridDimensions'][0] #2 + ord('x') - ord(axis)]
                proj = ds.proj(field,axis)
                frb = proj.to_frb(1,res)
                hdu = pyfits.PrimaryHDU(frb[field])
                hdulist = pyfits.HDUList([hdu])
                hdulist.writeto(outfile,clobber=True)
                print("wrote", outfile)
    def QUEB(self, frame, BoxSize=None):
        #ds = self.car.load(frame)
        frb_dir = "%s/%s/"%(self.directory,frbname)
        p49_QU2EB.QU2EB(frb_dir,frame,BoxSize=BoxSize)

    def EBslopes(self,frame, fit_range=None):
        EBSlopePower=p49_QU2EB.slopes_powers(frame,directory = self.directory + "/"+frbname, prefix=self.prefix, plot_format=self.plot_format,fit_range=fit_range )
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
        self.stuff['EBcycles'].append(frame)
        for ax in 'xyz':
            self.stuff['Eamp_%s'%ax  ].append(EBSlopePower['Eamp'][ax])
            self.stuff['Bamp_%s'%ax  ].append(EBSlopePower['Bamp'][ax])
            self.stuff['Eslope_%s'%ax].append(EBSlopePower['Eslope'][ax])
            self.stuff['Bslope_%s'%ax].append(EBSlopePower['Bslope'][ax])

    def read_queb(self,frame,ax='x'):
        def read_fits(fitname):
            d=None
            if os.path.exists(fitname):
                d=np.ascontiguousarray(pyfits.open(fitname)[0].data,dtype=np.double)
            return d


        frb_dir = "%s/%s"%(self.directory,frbname)
        product_dir = "%s/DD%04d.products"%(self.directory,frame)
        xd='DD'
        Df= "%s/%s%04d_density_%s.fits"%(frb_dir,xd,frame,ax)
        Hf= "%s/%s%04d_magnetic_field_strength_%s.fits"%(frb_dir,xd,frame,ax)
        Qf= "%s/%s%04d_Q%s.fits"%(frb_dir,xd,frame,ax)
        Uf= "%s/%s%04d_U%s.fits"%(frb_dir,xd,frame,ax)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,frame,ax)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,frame,ax)
        d=read_fits(Df)
        h=read_fits(Hf)
        q=read_fits(Qf)
        u=read_fits(Uf)
        e=read_fits(Ef)
        b=read_fits(Bf)
        ts=queb_snapshot(q,u,d,H=h,BoxSize=self.BoxSize,axis=ax,simulation=self,frame=frame)
        ts.compute_harmonic_products()
        return ts

    def image_fields(self,frame,axis='x',ts=None):
        if ts is None:
            ts=self.read_queb(frame,axis)
            ts.compute_harmonic_products()
        for name in [ 'T','H','Q','U','E','B']:
            fig,ax=plt.subplots(1,1)
            ax.clear()
            array = ts[name]
            proj=ax.imshow(array,origin='lower',interpolation='nearest')
            fig.colorbar(proj, ax=ax)
            outname = "%s_%04d_%s_%s.png"%(self.prefix,frame,name,axis)
            fig.savefig(outname)
            print(outname)
            plt.close(fig)
        return ts

    def make_spectra(self,frame):
        oober = st.short_oober(self.directory, frame=frame)
        st.MakeVelocitySpectra(oober,frame)
        st.MakeAccelSpectra(oober,frame)
        st.MakeVelocitySpectra(oober,frame)
        st.MakeMagneticSpectra(oober,frame)
        st.MakeDensitySpectra(oober,frame)

    def read_spectra(self,frame,fit_range=None,ax='x'):
        """reads spectra data from several sources.
        QUEB spectra are read from fits files, and the spectra are computed.
        Fluid quantity spectra are computed elsewhere (bug collins for the code)
        and stored in DD????.products.  
        These are plotted with plot, below.
        """
        frb_dir = "%s/%s"%(self.directory,frbname)
        xd='DD'
        Df= "%s/%s%04d_density_%s.fits"%(frb_dir,xd,frame,ax)
        Hf= "%s/%s%04d_magnetic_field_strength_%s.fits"%(frb_dir,xd,frame,ax)
        Qf= "%s/%s%04d_Q%s.fits"%(frb_dir,xd,frame,ax)
        Uf= "%s/%s%04d_U%s.fits"%(frb_dir,xd,frame,ax)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,frame,ax)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,frame,ax)
        d=np.ascontiguousarray(pyfits.open(Df)[0].data,dtype=np.double)
        h=np.ascontiguousarray(pyfits.open(Hf)[0].data,dtype=np.double)
        q=np.ascontiguousarray(pyfits.open(Qf)[0].data,dtype=np.double)
        u=np.ascontiguousarray(pyfits.open(Uf)[0].data,dtype=np.double)
        e=np.ascontiguousarray(pyfits.open(Ef)[0].data,dtype=np.double)
        b=np.ascontiguousarray(pyfits.open(Bf)[0].data,dtype=np.double)
        ts=queb_snapshot(q,u,d,H=h,BoxSize=self.BoxSize,axis=ax,simulation=self,frame=frame)
        ts.compute_harmonic_products()
        self.make_spectra(frame)
        ts.vspec=dt.dpy( "%s/DD%04d.products/power_velocity.h5"%(self.directory,frame) , ['k','power'])
        ts.aspec=dt.dpy( "%s/DD%04d.products/power_acceleration.h5"%(self.directory,frame) , ['k','power'])
        ts.dspec=dt.dpy( "%s/DD%04d.products/power_density.h5"%(self.directory,frame) , ['k','power'])
        ts.hspec=dt.dpy( "%s/DD%04d.products/power_magnetic.h5"%(self.directory,frame) , ['k','power'])
        ts.frame=frame
        ts.axis=ax
        ts.prefix=self.prefix
        #product_dir = "%s/DD%04d.products"%(self.directory,frame)
        #N_fft = dt.read_fft('%s/fft_density_%s.float32'%(product_dir,ax),'density')
        #ts['N_fft']=N_fft[:,:N_fft.shape[1]//2+1]
        return ts

    def plot_eb(self,ts,fname='TEST.png', **kwargs):
        """plot spectra extracted from e&b.
        the function 'dostuff' treats normalization and slicing."""

        def dostuff(arr):
            return arr[1:]#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])
        def dostuff2(arr):
            return arr[1:]/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])

        fig,ax=plt.subplots(1,1)
        k = ts.vspec[0]
        ylimits=dt.extents()
        xlim=dt.extents()
        pos_k=slice(None)
        this_k=(k/k[1])[pos_k]
        this_k=(k*2*np.pi)#[1:]
        this_k = 0.5*(this_k[1:]+this_k[:-1])
        ell=ts.lcent
        #this_ell=(ell/ell[1])[pos_k]
        this_ell=ell[:-1]#[1:] 

        rEB = ts['ClEB']/(ts['ClEE']*ts['ClBB'])**0.5
        lab=r'$r_{EB}=C_{\ell}^{EB}/\sqrt{C_{\ell}^{EE}C_{\ell}^{EE}}$'
        ax.plot( this_k,dostuff2(ts['aspec'][1]),c='k',marker='*', label=r'$a$');   ylimits(ts['aspec'][1][pos_k])# print('a lim',ylimits)
        ax.plot( this_ell,dostuff2(ts['ClEE']),        marker='*', label=r'$EE$',c='g'); ylimits(rEB)# print(ylimits)
        ax.plot( this_ell,dostuff2(ts['ClBB']),        marker='*', label=r'$BB$',c='b'); ylimits(rEB)# print(ylimits)
        ax.plot( this_ell,dostuff(rEB),                marker='*', label=r'$r_{EB}$',c='m'); ylimits(rEB)# print(ylimits)
        dt.axbonk(ax,xlabel='k/k1',ylabel=lab,xscale='log',yscale='log')
        #ax.set_xscale('symlog',linthreshx=1)
        #ax.set_xlim(xlim)
        #ax.set_ylim(1e-9,1e4)
        ax.set_yscale('symlog',linthreshy=1e-2)
        #ax.set_ylim(ylimits)
        ax.set_ylim(-1,1)
        #print("ylimits",ylimits)

        title="t2 %s n%04d %s"%(ts['prefix'],ts['frame'],ts['axis'])
        print("POOT",title)
        ax.legend(loc=0)
        ax.set_title(title)
        fig.savefig(fname)
        plt.close(fig)
        #ax.plot( this_k,dostuff(ts['aspec'][1][pos_k]),marker='*', label=r'$P(a)$');   ylimits(ts['aspec'][1][pos_k])# print('a lim',ylimits)
    def plot_many_spectra(self,ts,fname='TEST.png', **kwargs):
        def dostuff(arr):
            return arr[1:]#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])
        fig,ax=plt.subplots(1,1)
        k = ts['vspec'][0]
        ylimits=dt.extents()
        xlim=dt.extents()
        pos_k=slice(None)
        this_k=(k/k[1])[pos_k]
        this_k=(k*2*np.pi)#[1:]
        this_k = 0.5*(this_k[1:]+this_k[:-1])
        ax.plot( this_k,dostuff(ts['aspec'][1][pos_k]),marker='*', label=r'$P(a)$');   ylimits(ts['aspec'][1][pos_k])# print('a lim',ylimits)
        ax.plot( this_k,dostuff(ts['vspec'][1][pos_k]),marker='*', label=r'$P(v)$');   ylimits(ts['vspec'][1][pos_k])# print('v lim',ylimits)
        ax.plot( this_k,dostuff(ts['dspec'][1][pos_k]),marker='*', label=r'$P(\rho)$');ylimits(ts['dspec'][1][pos_k])# print('d lim',ylimits)
        ax.plot( this_k,dostuff(ts['hspec'][1][pos_k]),marker='*', label=r'$P(H)$');   ylimits(ts['hspec'][1][pos_k])# print('h lim',ylimits)
        ell=ts['lcent']
        #this_ell=(ell/ell[1])[pos_k]
        this_ell=ell[:-1]#[1:] 
        ax.plot( this_ell,dostuff(ts['ClEE'][pos_k]),marker='*', label=r'$C_{\ell}^{EE}$'); ylimits(ts['ClEE'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClBB'][pos_k]),marker='*', label=r'$C_{\ell}^{BB}$'); ylimits(ts['ClBB'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClTE'][pos_k]),marker='*', label=r'$ClTE$') ; ylimits(ts['ClTE'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClTT'][pos_k]),marker='*', label=r'$ClTT$') ; ylimits(ts['ClTT'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClHH'][pos_k]),marker='*', label=r'$ClHH$') ; ylimits(ts['ClHH'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClTB'][pos_k]),marker='*', label=r'$ClTB$') ; ylimits(ts['ClTB'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClEB'][pos_k]),marker='*', label=r'$ClEB$') ; ylimits(ts['ClEB'])# print(ylimits)
        dt.powerline(ax, this_ell[1]*4, this_ell[1]*10, 1, -5./3,c='k',label='-5/3')
        dt.powerline(ax, this_ell[1]*4, this_ell[1]*10, 1, -2.5,c='k',label='-2.5',linestyle='--')
        #ax.plot( this_ell, np.abs(ts['pork'])[pos_k],marker='*', label=r'pork') ; ylimits(ts['ClTE']); print(ylimits)
        #ts['ColumnDensity']= np.abs(cmbtools.harm2cl( ts['Eh'], ts['Deltal'],ts['lbins']))
        #ax.plot(this_ell,np.abs(ts['ColumnDensity'][pos_k]),c='k')
        #ax.plot( ell/ell[1], ts['ClTB'],marker='*', label=r'$ClEB$')
        xlim(this_k)
        xlim(this_ell)
        
        ax.legend(loc=1)
        ts.limits=ylimits


        #dt.axbonk(ax,xlabel='k/k1',ylabel='power',xscale='log',yscale='log')
        dt.axbonk(ax,xlabel='k/k1',ylabel='power',xscale='log',yscale='log')
        #ax.set_xscale('symlog',linthreshx=1)
        #ax.set_xlim(xlim)
        #ax.set_ylim(1e-9,1e4)
        ax.set_yscale('symlog',linthreshy=1e-28)
        #ax.set_ylim(ylimits)
        ax.set_ylim(-1,1)
        print("ylimits",ylimits)

        title="t2 %s n%04d %s"%(ts['prefix'],ts['frame'],ts['axis'])
        print("POOT",title)
        ax.set_title(title)
        fig.savefig(fname)
        plt.close(fig)
