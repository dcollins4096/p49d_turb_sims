from GL import *
import queb3
reload(queb3)
import torch
import simulation
import torch.nn.functional as F
from collections import defaultdict
import tqdm

def downsample_avg(x, M):
    if x.ndim == 2:   # [N, N]
        x = x.unsqueeze(0).unsqueeze(0)  # -> [1, 1, N, N]
        out = F.adaptive_avg_pool2d(x, (M, M))
        return out.squeeze(0).squeeze(0) # -> [M, M]
    elif x.ndim == 4: # [B, C, N, N]
        return F.adaptive_avg_pool2d(x, (M, M))
    else:
        raise ValueError("Input must be [N, N] or [B, C, N, N]")


def pull(simlist, size, N_per_frame, target_res = None,suffix="", rotate=False, los='xyz', half=None, fields='THQUEB'):
    output = []
    quan = defaultdict(list)
    adder=''
    if target_res:
        adder = "_down_%d"%target_res
    if rotate:
        adder += "_rot"
    whichhalf=''
    if half is not None:
        if half==0:
            whichhalf = '_first'
        if half==1:
            whichhalf='_second'
    oname = "p79d_subsets_S%d_N%d_%s_%s%s%s.h5"%(size,N_per_frame,los, adder, suffix,whichhalf)
    print(oname)
    for sim in simlist:
        print(sim)
        this_sim = simulation.corral[sim]
        this_sim.read_avg_quan()
        sl = slice(None)
        if half is not None:
            if half==0:
                sl = slice(0,len(this_sim.ann_frames)//2)
            if half==1:
                sl = slice(len(this_sim.ann_frames)//2, None)


        print(len(this_sim.ann_frames))
        for frame in this_sim.ann_frames[sl]:
            for ilos,this_los in enumerate(los):

                TEB={}
                for field in fields:
                    if field == 'T':
                        field_name = 'density_'
                    elif field == 'H':
                        field_name = 'magnetic_field_strength_'
                    elif field == 'P':
                        field_name = "H_POS_"
                    elif field == 'V':
                        field_name = 'velocity_centroid_'
                    elif field == 'S':
                        field_name = 'velocity_variance_'
                    else:
                        field_name = field
                    frb_name = "%s/DD%04d.products/DD%04d_%s%s.fits"%(this_sim.product_location,frame,frame,field_name,this_los)
                    arr = pyfits.open(frb_name)[0].data
                    TEB[field]=np.tile(arr,(2,2))
                    Nside = arr.shape[0]
                for n in range(N_per_frame):
                    corner = (np.random.random(2)*(Nside)).astype('int')
                    other = corner + np.array([size,size])
                    TEB1 = np.zeros([len(fields),size,size])
                    for nf,field in enumerate(fields):
                        TEB1[nf,:,:] = TEB[field][corner[0]:other[0], corner[1]:other[1]]
                    if rotate:
                        k=np.random.randint(0,4)
                        TEB1 = np.rot90(TEB1, k=k, axes=(1,2))
                    output.append(TEB1)

                    quan['frame'].append(frame)
                    iq = np.where( this_sim.quan_time['frames']==frame)[0][0]
                    quan['Ms_mean'].append( this_sim.Ms_mean)
                    quan['Ma_mean'].append( this_sim.Ma_mean)
                    quan['Ms_act'].append( this_sim.quan_time['vrms'][iq])
                    quan['Ma_act'].append( this_sim.quan_time['ma'][iq])
                    #print('ms mean %0.2f ms act %0.2f'%(this_sim.Ms_mean, quan['Ms_act'][-1]))
                    quan['los'].append(ilos)


    Nsubs = len(output)
    total = np.zeros([Nsubs, len(fields), size, size])
    total = torch.tensor(total,dtype=torch.float32)
    inds = np.arange(Nsubs)
    quan2 = {'Ms_mean':[],'Ma_mean':[],'Ms_act':[],'Ma_act':[],'los':[], 'frame':[]}
    print('randomize')
    import tqdm
    for n in tqdm.tqdm(inds):
        b = int((np.random.random()*len(output))//1)
        total[n,...] = torch.tensor(output.pop(b), dtype=torch.float32)
        for q in quan:
            quan2[q].append( quan[q].pop(b))
        #print(len(output))
    if target_res:
        print('downsample')
        total = downsample_avg(total, target_res)
    print('write')
    fptr = h5py.File(oname,'w')
    print(oname)
    fptr['subsets']=total
    for q in quan:
        fptr[q] = np.array(quan2[q])



    fptr.close()






