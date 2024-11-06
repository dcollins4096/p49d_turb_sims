#!/usr/bin/env python

import numpy as np
import sys
import os
import interpo
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s","--ms", dest="ms", action='store', help='sonic mach', type="float")
parser.add_option("-a","--ma", dest="ma", action='store', help='alfven mach', type="float")
parser.add_option("-n","--name", dest="name", action='store', help='job name', type="string")
parser.add_option("-d","--dimension", dest="dims", action='store', help='top grid side length', type="int")

density=1
N_dynamical_times = 2
options, args = parser.parse_args()

mach_nom = options.ms
alfmach = options.ma

mach = interpo.get_ms_nom(mach_nom,alfmach)





if np.abs(alfmach) < 1e-4:
    B=0
else:
    B = mach*np.sqrt(density)/alfmach*np.sqrt(4*np.pi)
#mach1d = mach/np.sqrt(3)
mach1d = mach
args = {}
args['TopGridDimension'] = options.dims
args['mach1d'] = mach1d
args['Bfield']=B
tdyn = 0.5/mach
args['tstop']=N_dynamical_times*tdyn
args['dt'] = tdyn/10


outdir = "%s_Ms%0.1f_Ma%0.1f_%d"%(options.name,mach,alfmach,options.dims)

if not os.path.exists(outdir):
    os.makedirs(outdir)

fname = outdir + "/" + options.name +".enzo"
print(args)

import jinja2
loader=jinja2.FileSystemLoader('.')
env = jinja2.Environment(loader=loader)
template = env.get_template('Template.enzo')
foutptr = open(fname,'w')
foutptr.write( template.render(**args))
foutptr.close()

bargs={}
bargs['jobname']=options.name
bargs['enzoname']=options.name+".enzo"
bfname = outdir + "/Job.sbatch"

loader=jinja2.FileSystemLoader('.')
env = jinja2.Environment(loader=loader)
template = env.get_template('Template.sbatch')
foutptr = open(bfname,'w')
foutptr.write( template.render(**bargs))
foutptr.close()

import shutil

shutil.copy2("enzo.exe", "%s/enzo.exe"%outdir)


#end
