from GL import *
import sim_colors
import queb3
import yt
reload(queb3)
#Note that the P49:Mach<5 runs were done with an inline version of the averaging tool.

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
    def ytcall(self, all_data):
        #pdb.set_trace()
        Q = all_data[self.quantity]
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
    def ytcall(self,all_data):
        d = all_data['density']
        Bi = all_data[self.quantity]
        Q = Bi/np.sqrt(d)
        self.avg +=  Q.mean()
        self.var += (Q**2).mean()
        self.N1+=1
        self.N2+=1
    def __call__(self,grid):
        d = grid['Density']
        Bi = grid[self.quantity]
        Q = Bi/np.sqrt(d)
        self.avg +=  Q.mean()
        self.var += (Q**2).mean()
        self.N1+=1
        self.N2+=1

import re
def parse_athena_meta(fname):
    fptr = open(fname)
        #<Time Value="0.412533"/> 
    rrr = re.compile(r"\s*<Time Value=\"(.*)\"/>")
    lines = fptr.readlines()
    fptr.close()
    for line in lines:
        match = rrr.match(line)
        if match:
            time = float(match.group(1))
            break
    return time


def make_quan_athena(directory,frame, out_directory=None, clobber=False, sim='SIM', do_magnetic=True):
    #for athena.

    input_file = "%s/parthenon.prim.%05d.phdf"%(directory, frame)
    input_xml = "%s/parthenon.prim.%05d.phdf.xdmf"%(directory, frame)
    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    outname_flag = "%s/DD%04d.products/data%04d.AverageQuantities.h5.flag"%(out_directory,frame,frame)
    #print(outname)
    if os.path.exists(outname) and clobber==False:
        print("File exists, skipping", outname)
        return 0
    if os.path.exists(outname_flag):
        #being worked on.
        return 0

    time = nar([parse_athena_meta(input_xml)])
    print("Quan on frame",frame)
    submarine={}
    submarine['density']=meanie('density','density')
    submarine['vx']=meanie('vx','velocity_x')
    submarine['vy']=meanie('vy','velocity_y')
    submarine['vz']=meanie('vz','velocity_z')
    if do_magnetic:
        submarine['bx']=meanie('bx','magnetic_field_x')
        submarine['by']=meanie('by','magnetic_field_y')
        submarine['bz']=meanie('bz','magnetic_field_z')
        submarine['alf_x']=alfv_meanie('alf_x','magnetic_field_x')
        submarine['alf_y']=alfv_meanie('alf_y','magnetic_field_y')
        submarine['alf_z']=alfv_meanie('alf_z','magnetic_field_z')

    ds = yt.load(input_file)
    all_data = ds.all_data()
    for sub in submarine:
        submarine[sub].ytcall(all_data)
    #do all averages

    for sub in submarine:
        submarine[sub].finish()


    parent_dir = os.path.dirname(outname)
    if not os.path.exists(parent_dir):
        grandparent_dir = os.path.dirname(parent_dir)
        if not os.path.exists(grandparent_dir):
            os.mkdir(grandparent_dir)
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


    
def make_quan(directory,frame, out_directory=None, clobber=False, sim='SIM', do_magnetic=True):
    #for enzo.

    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    outname_flag = "%s/DD%04d.products/data%04d.AverageQuantities.h5.flag"%(out_directory,frame,frame)
    outname_short = "./%s/DD%04d.products/data%04d.AverageQuantities.h5"%(sim,frame,frame)
    parent_dir = os.path.dirname(outname)
    if not os.path.exists(parent_dir):
        grandparent_dir = os.path.dirname(parent_dir)
        if not os.path.exists(grandparent_dir):
            os.mkdir(grandparent_dir)
        os.mkdir(parent_dir)
    #print(outname)
    if (queb3.check_finished(outname_short) or os.path.exists(outname) or os.path.exists(outname_flag)) and clobber==False:
        print("File exists, skipping", outname)
        return 0
    print("Quan on frame",frame)
    fptr = open(outname_flag,'w')
    fptr.close()
    submarine={}
    submarine['density']=meanie('density','Density')
    submarine['vx']=meanie('vx','x-velocity')
    submarine['vy']=meanie('vy','y-velocity')
    submarine['vz']=meanie('vz','z-velocity')
    if do_magnetic:
        submarine['bx']=meanie('bx','Bx')
        submarine['by']=meanie('by','By')
        submarine['bz']=meanie('bz','Bz')
        submarine['alf_x']=alfv_meanie('alf_x','Bx')
        submarine['alf_y']=alfv_meanie('alf_y','By')
        submarine['alf_z']=alfv_meanie('alf_z','Bz')


    file_glob = "%s/DD%04d/data%04d.cpu*"%(directory,frame,frame)
    file_list=sorted(glob.glob(file_glob))

    #do all averages
    total=len(file_list)
    for n,fname in enumerate(file_list):
        #print("     ",fname, "%d/%d"%(n,total))
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


    










