
import xml.etree.ElementTree as ET
import compute_avg_quantities as caq
reload(caq)



#x = ET.parse('/anvil/scratch/x-ux454321/p83_turbulence/Athena/aa_Ms2.0_Ma0.5_512/parthenon.prim.00033.phdf.xdmf')
x = caq.parse_athena_meta('/anvil/scratch/x-ux454321/p83_turbulence/Athena/aa_Ms2.0_Ma0.5_512/parthenon.prim.00033.phdf.xdmf')
