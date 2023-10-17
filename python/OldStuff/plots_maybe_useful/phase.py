from GL import *
import yt
import sim_colors

def add_v2thing(obj):
    def v2(field,data):
        B0 = data.get_field_parameter('B0')
        return 1.5*data['velocity_magnitude']/B0
    obj.add_field("V2",v2,validators=[yt.ValidateParameter('B0')],sampling_type='cell',units='cm/(G*s)')
    def bfluct(field,data):
        B = data.get_field_parameter('B')
        if 0==B:
            return data['magnetic_field_strength']
        bx = data['magnetic_field_x']-B[0]
        by = data['magnetic_field_y']-B[1]
        bz = data['magnetic_field_z']-B[2]
        b = np.sqrt( bx**2 + by**2 + bz**2 )
        return b
    obj.add_field("bfluct",bfluct,validators=[yt.ValidateParameter('B')],sampling_type='cell', units='Gauss')


for sim in sim_colors.simlist[:1]:
    for frame in sim_colors.frames[sim][0:1]:
        fname = "%s/%s/DD%04d/data%04d"%(sim_colors.cloudbreak_base, sim, frame, frame)
        print(os.path.exists(fname))
        ds = yt.load(fname)
        add_v2thing(ds)
        #reg = ds.region([0.05]*3,[0.0]*3,[0.1]*3)
        reg=ds.all_data()
        bx = reg['magnetic_field_x']
        by = reg['magnetic_field_y']
        bz = reg['magnetic_field_z']
        B0x = bx.mean()
        B0y = by.mean()
        B0z = bz.mean()
        B0 = np.sqrt( B0x**2 + B0y**2 + B0z**2 )
        reg.set_field_parameter('B0',B0)
        reg.set_field_parameter('B', [B0x, B0y, B0z])
        phase = yt.create_profile(reg,['bfluct','V2'],'cell_volume', weight_field=None)
        pp = yt.PhasePlot.from_profile(phase)
        pp.set_xlabel('bfluct')
        pp.set_ylabel('v^2/B')
        pp.save("%s/phase_bfluct_V2_%s_n%04d"%(plotdir, sim,frame))




