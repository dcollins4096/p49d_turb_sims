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
class Ekin_meanie(meanie):
    def __init__(self,name,enzo_quantity):
        super().__init__(name,enzo_quantity)
    def __call__(self,grid):
        d = grid['Density'][()]
        #dx = grid['x-acceleration'][()]
        #dy = grid['y-acceleration'][()]
        #dz = grid['z-acceleration'][()]
        vx = grid['x-velocity'][()]
        vy = grid['y-velocity'][()]
        vz = grid['z-velocity'][()]
        Q = (0.5*d*(vx**2+vy**2+vz**2))
        self.avg +=  Q.sum()
        self.var += (Q**2).sum()
        self.N1+=Q.size
        self.N2+=Q.size
class Edot_meanie(meanie):
    def __init__(self,name,enzo_quantity):
        super().__init__(name,enzo_quantity)
    def __call__(self,grid):
        d = grid['Density'][()]
        dx = grid['x-acceleration'][()]
        dy = grid['y-acceleration'][()]
        dz = grid['z-acceleration'][()]
        vx = grid['x-velocity'][()]
        vy = grid['y-velocity'][()]
        vz = grid['z-velocity'][()]
        Q = (0.5*d*(vx*dx+vy*dy+vz*dz))
        self.avg +=  Q.sum()
        self.var += (Q**2).sum()
        self.N1+=Q.size
        self.N2+=Q.size

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

def bulk_viscosity_estimate(directory,frame,out_directory=None,sim='SIM', clobber=False):
    #outname = "%s/DD%04d.products/data%04d.BulkViscosity.h5"%(out_directory,frame,frame)
    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    #print("Bulk on frame",frame)
    print("Add bulk to",outname)
    optr = h5py.File(outname, 'r+')
    if 'vorticity_avg' in optr and not clobber:
        print( "Exists.  Skipping")
        optr.close()
        return
    optr.close()
    ds_name = "%s/DD%04d/data%04d"%(directory,frame,frame)
    ds = yt.load(ds_name)
    #yt.ProjectionPlot(ds,0,'vorticity_magnitude').save('%s/omega'%plot_dir)
    ad = ds.all_data()
    #ad = ds.region([0.5,0.5,0.5],[0.25,0.25,0.25],[0.75,0.75,0.75])
    omega2 = ad['vorticity_magnitude']**2
    div2   = ad['velocity_divergence']**2
    optr = h5py.File(outname, 'r+')
    try:
        print('means')
        optr['vorticity_avg'] = nar([omega2.mean().v])
        optr['vorticity_std'] = nar([omega2.std().v])
        optr['divergence_avg'] =nar([div2.mean().v])
        optr['divergence_std'] =nar([div2.std().v])
    except:
        raise
    finally:
        optr.close()

def make_ekin(directory,frame,out_directory=None,sim='SIM', clobber=False):
    #outname = "%s/DD%04d.products/data%04d.BulkViscosity.h5"%(out_directory,frame,frame)
    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    #print("Bulk on frame",frame)
    print("Add Ekin to",outname)
    optr = h5py.File(outname, 'r+')
    if 'Ekin' in optr and not clobber:
        print( "Exists.  Skipping")
        optr.close()
        return
    optr.close()
    ds_name = "%s/DD%04d/data%04d"%(directory,frame,frame)
    ds = yt.load(ds_name)
    #yt.ProjectionPlot(ds,0,'vorticity_magnitude').save('%s/omega'%plot_dir)
    ad = ds.all_data()
    #ad = ds.region([0.5,0.5,0.5],[0.25,0.25,0.25],[0.75,0.75,0.75])

    print('read vel')
    vx = ad['x-velocity'].v
    print('read vel')
    vy = ad['y-velocity'].v
    print('read vel')
    vz = ad['z-velocity'].v
    print('read den')
    rho = ad['density'].v

    Ekinetic = (0.5*rho*(vx**2+vy**2+vz**2)).sum()/rho.size

    optr = h5py.File(outname, 'r+')
    try:
        optr['Ekin'] = nar([Ekinetic]) 
    except:
        raise
    finally:
        optr.close()

def make_edot(directory,frame,out_directory=None,sim='SIM', clobber=False):
    #outname = "%s/DD%04d.products/data%04d.BulkViscosity.h5"%(out_directory,frame,frame)
    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    #print("Bulk on frame",frame)
    print("Add Edot to",outname)
    optr = h5py.File(outname, 'r+')
    if 'Edot' in optr and not clobber:
        del optr['Edot']
        optr.close()
        #return
    optr.close()
    ds_name = "%s/DD%04d/data%04d"%(directory,frame,frame)
    ds = yt.load(ds_name)
    #yt.ProjectionPlot(ds,0,'vorticity_magnitude').save('%s/omega'%plot_dir)
    ad = ds.all_data()
    #ad = ds.region([0.5,0.5,0.5],[0.25,0.25,0.25],[0.75,0.75,0.75])
    print('read drive')
    dx = ad['x-acceleration'].v
    print('read drive')
    dy = ad['y-acceleration'].v
    print('read drive')
    dz = ad['z-acceleration'].v

    print('read vel')
    vx = ad['x-velocity'].v
    print('read vel')
    vy = ad['y-velocity'].v
    print('read vel')
    vz = ad['z-velocity'].v
    
    print('read den')
    rho = ad['density'].v
    eta = ds['DrivingEfficiency']

    Edot = (rho*(dx*vx+dy*vy+dz*vz)*eta).sum()/rho.size
    Ekinetic = (0.5*rho*(vx**2+vy**2+vz**2)).sum()/rho.size
    Driving = (0.5*rho*(dx**2+dy**2+dz**2)).sum()/rho.size
    #pdb.set_trace()

    optr = h5py.File(outname, 'r+')
    try:
        if 'Edot' not in optr:
            optr['Edot'] = nar([Edot]) 
        if 'Ekin' not in optr:
            optr['Ekin'] = nar([Ekinetic]) 
    except:
        raise
    finally:
        optr.close()



def make_edot_faster(directory,frame,out_directory=None,sim='SIM', clobber=False):
    #outname = "%s/DD%04d.products/data%04d.BulkViscosity.h5"%(out_directory,frame,frame)
    outname = "%s/DD%04d.products/data%04d.AverageQuantities.h5"%(out_directory,frame,frame)
    #print("Bulk on frame",frame)
    print("Add Edot to",outname)
    optr = h5py.File(outname, 'r+')
    if 'Edot_avg' in optr:
        print("Exists.  Skip.")
        optr.close()
        return
        
    #for Q in ['Ekin','Edot_avg','Edot_std','Ekin_avg','Ekin_std']:
    #    if Q in optr and not clobber:
    #        del optr[Q]
    #        #return
    optr.close()
    submarine = {}
    submarine['Ekin'] = Ekin_meanie('Ekin','Ekin')
    submarine['Edot'] = Edot_meanie('Edot','Edot')

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

    optr = h5py.File(outname, 'r+')
    try:
        for sub in submarine:
            name="%s_avg"%submarine[sub].name
            if name not in optr:
                optr[name]=nar([submarine[sub].avg])
            name="%s_std"%submarine[sub].name
            if name not in optr:
                optr[name]=nar([submarine[sub].std])
    except:
        raise
    finally:
        optr.close()



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


    










