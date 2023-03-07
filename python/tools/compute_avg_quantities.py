from GL import *
import sim_colors
#Note that the Mach<5 runs were done with an inline version of the averaging tool.

class meanie():
    def __init__(self,name,enzo_quantity):
        self.name=name
        self.quantity=enzo_quantity
        self.N1=0
        self.N2=0
        self.avg = 0
        self.var = 0
        self.std = 0
    def __call__(self, grid):
        #pdb.set_trace()
        Q = grid[self.quantity][()]
        self.avg +=  Q.mean()
        self.var += (Q**2).mean()
        self.N1+=1
        self.N2+=1
    def finish(self):
        self.avg/=self.N1
        self.var /= self.N2
        self.std = np.sqrt( self.var - self.avg**2)

class alfv_meanie(meanie):
    def __init__(self,name,enzo_quantity):
        super().__init__(name,enzo_quantity)
    def __call__(self,grid):
        d = grid['Density']
        Bi = grid[self.quantity]
        Q = Bi/np.sqrt(d)
        self.avg +=  Q.mean()
        self.var += (Q**2).mean()
        self.N1+=1
        self.N2+=1

def make_quan(directory,frame, out_directory=None, clobber=False):
    if out_directory is None:
        out_directory=directory

    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    #print(outname)
    if os.path.exists(outname) and clobber==False:
        print("File exists, skipping", outname)
        return 0
    print("Quan on frame",frame)
    submarine={}
    submarine['bx']=meanie('bx','Bx')
    submarine['by']=meanie('by','By')
    submarine['bz']=meanie('bz','Bz')
    submarine['vx']=meanie('vx','x-velocity')
    submarine['vy']=meanie('vy','y-velocity')
    submarine['vz']=meanie('vz','z-velocity')
    submarine['alf_x']=alfv_meanie('alf_x','Bx')
    submarine['alf_y']=alfv_meanie('alf_y','By')
    submarine['alf_z']=alfv_meanie('alf_z','Bz')
    submarine['density']=meanie('density','Density')


    file_glob = "%s/DD%04d/data%04d.cpu*"%(directory,frame,frame)
    file_list=sorted(glob.glob(file_glob))

    #do all averages
    total=len(file_list)
    for n,fname in enumerate(file_list):
        print("     ",fname, "%d/%d"%(n,total))
        fptr = h5py.File(fname,'r')
        try:
            for grid in fptr:
                if grid.startswith('Meta'):
                    continue
                #pdb.set_trace()
                for sub in submarine:
                    submarine[sub](fptr[grid])

        except:
            raise
        finally:
            fptr.close()

    for sub in submarine:
        submarine[sub].finish()

    param_name = "%s/DD%04d/data%04d"%(directory,frame,frame)
    pptr=open(param_name,'r')
    for line in pptr.readlines():
        if line.startswith('InitialTime'):
            spl = line.split('=')
            time = nar([float(spl[1])])
            break
    pptr.close()




    outname = "%s/DD%04d/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    parent_dir = os.path.dirname(outname)
    if not os.path.exists(parent_dir):
        os.mkdir(parent_dir)

    optr = h5py.File(outname,'w')
    try:
        optr['time'] = time
        for sub in submarine:
            optr["%s_avg"%submarine[sub].name]=nar([submarine[sub].avg])
            optr["%s_std"%submarine[sub].name]=nar([submarine[sub].std])
    except:
        raise
    finally:
        optr.close()


    










