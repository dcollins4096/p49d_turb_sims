from GL import *
from collections import defaultdict

def return_average_quantities(file_list=[], all_quan = None):
    if all_quan is None:
        #all_quan = defaultdict(lambda: list())
        all_quan = defaultdict(list)
    else:
        all_quan = list(all_quan)
    for this_name in file_list:
        if not os.path.exists(this_name):
            print("No such file, skipping: %s"%this_name)
            continue
        fptr = h5py.File(this_name,'r')
        for key in fptr:
            mything =  list(fptr[key][:].flatten())
            all_quan[key] += mything
        try:
            pass
        except:
            raise
        finally:
            fptr.close()
    for key in all_quan:
        all_quan[key] = np.array( all_quan[key]).flatten()
    return all_quan

def files_from_output(path_to_output_log, suffix=".AverageQuantities.h5"):
    if 'OutputLog' not in path_to_output_log.split("/"):
        output_log = path_to_output_log + "/OutputLog"
        path = path_to_output_log
    else:
        output_log=path_to_output_log
        path = "/".join(output_log.split("/")[:-1])

    fptr = open(output_log,'r')
    lines = fptr.readlines()
    fptr.close()
    files = []
    for line in lines:
        files.append( path+line.split()[2][1:]+suffix)
    return files

def all_quan_from_outputlog(output_log):
    file_list= files_from_output(output_log)
    quan = return_average_quantities(file_list)
    return quan

def all_quan_from_taxi(car):
    dumb=[]
    file_list=[]
    for n in car.frame_dict:
        file_list.append("%s/%s.AverageQuantities.h5"%(car.directory,car.frame_dict[n]['dsname']))
    all_quan = return_average_quantities(file_list)
    return all_quan

def get_quantities_and_rms(output_log):
    quan = all_quan_from_outputlog(output_log)
    
    mytime=quan['time']
    mytime -= mytime[0]
    mytime /= mytime[-1]
    quan['mytime']=mytime

    vx2 = quan['vx_std']**2
    vy2 = quan['vy_std']**2
    vz2 = quan['vz_std']**2
    quan['vrms'] = np.sqrt(vx2 + vy2 + vz2)
    quan['Ms'] = quan['vrms']

    vax2 = quan['alf_x_std']**2
    vay2 = quan['alf_y_std']**2
    vaz2 = quan['alf_z_std']**2
    quan['varms']= np.sqrt(vax2 + vay2 + vaz2)

    Bx = quan['bx_avg']
    By = quan['by_avg']
    Bz = quan['bz_avg']
    Btot = np.sqrt( Bx**2+By**2+Bz**2)
    quan['Btot'] = Btot

    quan['Ma'] = quan['vrms']/quan['Btot']**2

    bx = quan['bx_std']**2
    by = quan['by_std']**2
    bz = quan['bz_std']**2
    brms = np.sqrt( bx + by + bz)
    #ax_b2v2.scatter( vrms.mean(), 1/(vrms.mean()/brms.mean()/Btot.mean()), c=sim_colors.color[sim],marker=sim_colors.marker[sim]) 
    quan['brms'] = brms
    quan['also_brms'] = quan['vrms']**2/quan['Btot']
    return quan

