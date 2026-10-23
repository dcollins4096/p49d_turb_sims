
from GL import *
import gc

import simulation

def ppv_cube(density, velocity, axis, vbins):
    """Density-weighted position-position-velocity cube along one axis.
    ppv[i,j,k] = (1/N_los) * sum of density over LOS cells with velocity in channel k,
    so ppv.sum(axis=-1) is the LOS mean density (same normalization as make_frbs),
    minus any mass outside the velocity range.
    The two image axes are the remaining spatial axes in xyz order."""
    axnum = 'xyz'.index(axis)
    nv = len(vbins)-1
    #put the line of sight last, so the image axes come first.
    rho = np.moveaxis(density, axnum, -1)
    vel = np.moveaxis(velocity, axnum, -1)
    n0, n1, nlos = rho.shape
    chan = np.digitize(vel, vbins) - 1
    ok = (chan >= 0) & (chan < nv)
    pix = np.broadcast_to( (np.arange(n0)[:,None]*n1 + np.arange(n1)[None,:])[:,:,None], rho.shape)
    index = pix[ok]*nv + chan[ok]
    ppv = np.bincount(index, weights=rho[ok], minlength=n0*n1*nv).reshape(n0,n1,nv)/nlos
    clipped = rho[~ok].sum()/rho.sum()
    return ppv, clipped

def make_ppv(simlist, axes='xyz', vmax=None, nsigma=5.0, nv=128, frames=None):
    """Write DDNNNN_ppv_<axis>.h5 to product_location/DDNNNN.products/ for each frame.
    nv velocity channels span (-vmax,vmax) in code units.  If vmax is None, it is set per frame
    and axis to |mean| + nsigma*rms of that velocity component, so channels track Mach number;
    the bins are saved with each cube.  No thermal broadening is applied."""
    for sim in simlist:
        this_sim = simulation.corral[sim]
        these_frames = this_sim.ann_frames if frames is None else frames
        for frame in these_frames:
            frame_dir = "%s/DD%04d.products"%(this_sim.product_location, frame)
            outnames = {axis:"%s/DD%04d_ppv_%s.h5"%(frame_dir, frame, axis) for axis in axes}
            todo = [axis for axis in axes if not os.path.exists(outnames[axis])]
            if len(todo) == 0:
                print("exists, continue:", list(outnames.values()))
                continue
            os.makedirs(frame_dir, exist_ok=True)
            ds = this_sim.load_ds(frame)
            cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
            density = cg['density'].v
            for axis in todo:
                velocity = cg['velocity_%s'%axis].v
                this_vmax = np.abs(velocity.mean()) + nsigma*velocity.std() if vmax is None else vmax
                vbins = np.linspace(-this_vmax, this_vmax, nv+1)
                vcen = 0.5*(vbins[:-1]+vbins[1:])
                ppv, clipped = ppv_cube(density, velocity, axis, vbins)
                del velocity
                cg.field_data.pop(('gas','velocity_%s'%axis), None)
                if clipped > 1e-3:
                    print("Warning: %0.2e of the mass is outside +-%0.2f"%(clipped, this_vmax))
                #write to a temp name and rename, so a crash never leaves a partial file.
                tmpname = outnames[axis] + ".tmp"
                h5ptr = h5py.File(tmpname, 'w')
                h5ptr.create_dataset('ppv', data=ppv.astype('float32'), compression='gzip')
                h5ptr['vbins'] = vbins
                h5ptr['vcen'] = vcen
                h5ptr['ppv'].attrs['axis'] = axis
                h5ptr['ppv'].attrs['image_axes'] = ''.join(a for a in 'xyz' if a != axis)
                h5ptr['ppv'].attrs['norm'] = 'sum over channels = LOS mean density'
                h5ptr['ppv'].attrs['clipped_mass_fraction'] = clipped
                h5ptr['ppv'].attrs['sim'] = sim
                h5ptr['ppv'].attrs['frame'] = frame
                h5ptr.close()
                os.replace(tmpname, outnames[axis])
                print("wrote", outnames[axis])
            del density, cg
            gc.collect()
