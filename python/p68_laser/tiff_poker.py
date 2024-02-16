
from dtools.starter1 import *

from tifffile import imread

base_dir="/Users/dcollins/Dropbox/RESEARCH5/Paper68/Data_analysis/Raw_radiographs/play"
i1="TD_TC090-124_HGXD_IMAGE_N220712-002-999_DROOP_CORR_422421128478532_20220907115837893.tif"
i2="TD_TC090-124_HGXD_IMAGE_N220713-001-999_DROOP_CORR_745405306609760_20220907115905225.tif"
i3="TD_TC090-124_HGXD_IMAGE_N220714-001-999_DROOP_CORR_766909572834217_20220907115956892.tif"
fnames=[i1,i2,i3]

plotdir="plots_to_sort"

class viewer():
    def __init__(self,which=0):
        self.fname = "%s/%s"%(base_dir,fnames[which])
        self.all_data = imread(self.fname)
        self.XX = np.arange(self.all_data.shape[1])
        self.YY = np.arange(self.all_data.shape[0])
        self.X1, self.Y1 = np.meshgrid(self.XX,self.YY)

        print(self.all_data.shape)
        print(self.X1.shape)

    def image1(self):

        self.guess_scale()
        vmin,vmax=self.minmax
        #mpl.image.imsave('plots_to_sort/full_uniform.png', self.all_data, cmap='viridis', vmin=self.minmax[0],vmax=self.minmax[1])
        norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        fig,axes=plt.subplots(1,1)
        #axes.imshow(self.all_data,norm=norm)
        axes.pcolormesh(self.X1,self.Y1,self.all_data,norm=norm)
        axes.set_aspect('equal')
        fig.savefig('plots_to_sort/image0')




    def guess_scale(self):
        ad=self.all_data
        test = ad[ad.shape[0]//4:int(3/4*ad.shape[0]),ad.shape[1]//4:int(3/4*ad.shape[1])]
        self.minmax=np.array([test[test>0].min(),test.max()])



    def image(self,a=None,b=None,c=None,d=None,vmin=None,vmax=None):
        fig,axes=plt.subplots(2,2,figsize=(12,12))
        #ax0=axes
        ax0=axes[0][0];ax1=axes[0][1]
        ax2=axes[1][0];ax3=axes[1][1]
        for ax in axes.flatten():
            ax.set_aspect('equal')
        norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        S1 = slice(c,d)
        S2 = slice(a,b)
        ax0.pcolormesh(self.X1,self.Y1,self.all_data,norm=norm)
        ax0.plot([a,b,b,a,a],[c,c,d,d,c],c='r')
        TheX,TheY,TheZ=self.X1[S1,S2],self.Y1[S1,S2],self.all_data[S1,S2]
        ax1.pcolormesh(TheX,TheY,TheX,norm=norm)

        for X in np.arange(a,a+100,5):
            ax0.axvline(X)
            ax1.axvline(X)
            TheX_s,TheY_s= TheY[X-a,:], TheZ[X-a,:]
            ax2.plot(TheX_s,TheY_s)



        fig.savefig("%s/image"%(plotdir))






