#
# Import NumPy for array handling
#
from GL import *
import numpy as np
import math
import yt
import numpy as np
import sim_colors 
from yt.extensions.astro_analysis.radmc3d_export.api import RadMC3DWriter, RadMC3DSource
reload(sim_colors)
import A2.write_rmc_fields as write_rmc_fields
reload(write_rmc_fields)

class rmc_filemaker():
    def __init__(self, repo_dir=None,output_dir=None, dust_model = 'simple_wrong'):
        self.repo_dir=repo_dir
        self.output_dir=output_dir
        self.dust_model = dust_model
        self.small_data_dir = "%s/data_files_small"%self.repo_dir
    def write_dust_parameters(self):
        if self.dust_model == 'simple_wrong':
            self.write_wrong_dust()
        else:
            print("Dust model not defined", self.dust_model)
    def write_wrong_dust(self):

        #
        # Make the wavelength_micron.inp file
        #
        lam1     = 0.1e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size
        #
        # parse the dustkapscatmat file to check how many wavelengths it has
        # (necessary for creating the mock alignment factor file)
        #
        with open('%s/dustkapscatmat_pyrmg70.inp'%self.small_data_dir,'r') as f:
            for _ in range(7): f.readline()
            dustnf = int(f.readline())
            f.readline()
            dustfreq = np.zeros(dustnf)
            f.readline()
            for inu in range(dustnf):
                s=f.readline().split()
                dustfreq[inu] = float(s[0])
        #
        # Now make a mock alignment factor model. This is ONLY FOR TESTING.
        # 
        nrang = 20
        muang = np.linspace(1.e0,0.e0,nrang)
        eta   = np.arccos(muang)*180./math.pi
        orth  = np.zeros(nrang) + 1.e0
        ampl  = 0.5
        para  = ( 1.e0 - ampl*np.cos(muang*math.pi) ) / ( 1.e0 + ampl)
        #
        # Write the wavelength file
        #
        with open('%s/wavelength_micron.inp'%(self.output_dir),'w+') as f:
            f.write('%d\n'%(nlam))
            for value in lam:
                f.write('%13.6e\n'%(value))
        #
        # Write the stars.inp file
        #
        if 0:
            with open('stars.inp','w+') as f:
                "We don't need any stars, we'll set the temperature by hand"
                pass

        #
        # Dust opacity control file
        #
        with open('%s/dustopac.inp'%(self.output_dir),'w+') as f:
            f.write('2               Format number of this file\n')
            f.write('1               Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('20              Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('pyrmg70         Extension of name of dustkappa_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')

        #
        # Dust alignment data
        #
        with open('%s/dustkapalignfact_pyrmg70.inp'%(self.output_dir),'w+') as f:
            f.write('1\n')
            f.write('%d\n'%(dustnf))
            f.write('%d\n\n'%(nrang))
            for value in dustfreq:
                f.write('%13.6e\n'%(value))
            f.write('\n')
            for value in eta:
                f.write('%13.6e\n'%(value))
            f.write('\n')
            for inu in range(dustnf):
                for imu in range(nrang):
                    f.write('%13.6e %13.6e\n'%(orth[imu],para[imu]))
                f.write('\n')

        #
        # Write the radmc3d.inp control file
        #
        nphot    = 1000000
        with open('%s/radmc3d.inp'%(self.output_dir),'w+') as f:
            f.write('nphot = %d\n'%(nphot))
            f.write('scattering_mode_max = 4\n')
            f.write('alignment_mode = -1\n')
            f.write('iranfreqmode = 1\n')
    def write_radmc_from_enzo(self, enzo_fname, output_dir=None):

        #
        # Write the grid file
        #

        ds = yt.load(enzo_fname)


        writer = RadMC3DWriter(ds)
        writer.write_amr_grid()
        os.mv


        #
        # Write the density file
        #
        dust_to_gas = 0.01
        def _DustDensity(field, data):
            return dust_to_gas * data["density"]
        ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3",sampling_type='cell')

        def _DumbTemp(field, data):
            return ds.arr(np.ones_like(data["density"]), 'K')*10
        ds.add_field("DumbTemp", function=_DumbTemp, units="K", sampling_type='cell')


        #radmc3d complains if the alignment vector is >1, even by 1e-16.
        #eps_val keeps it from complaining.
        eps_val = 1e-16
        def _bx_hat(field, data):
            eps = data.ds.quan(eps_val,'gauss')
            output = (data['magnetic_field_x']/(data['magnetic_field_strength']+eps))
            return output
        ds.add_field("bx_hat", function=_bx_hat, units="dimensionless",sampling_type='cell')
        def _by_hat(field, data):
            eps = data.ds.quan(eps_val,'gauss')
            output = data['magnetic_field_y']/(data['magnetic_field_strength']+eps)
            return output
        ds.add_field("by_hat", function=_by_hat, units="dimensionless",sampling_type='cell')
        def _bz_hat(field, data):
            eps = data.ds.quan(eps_val,'gauss')
            output = data['magnetic_field_z']/(data['magnetic_field_strength']+eps)
            return output
        ds.add_field("bz_hat", function=_bz_hat, units="dimensionless",sampling_type='cell')
        write_rmc_fields.write_rmc_fields(writer,("gas","dust_density"), "dust_density.inp")
        write_rmc_fields.write_rmc_fields(writer,"DumbTemp", "dust_temperature.dat")
        write_rmc_fields.write_rmc_fields(writer,["bx_hat","by_hat","bz_hat"], "grainalign_dir.inp")
