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
