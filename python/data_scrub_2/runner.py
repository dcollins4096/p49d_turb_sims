from  GL import *

import simulation
import simulation_info.all_sims

import data_scrub_2.make_quan as caq
reload(caq)
import data_scrub_2.make_2d_spectra as m2d
import data_scrub_2.make_3d_spectra as m3d
import data_scrub_2.make_all_frb as maf
#import data_scrub_2.image_frb as imf
import plots.P3_all_spectra as p3
import plots.P1_plot_quan as p1
import p49_plots.rework_spectra as rework
#this_simname="1_1"
#this_simname="4s_dev"
#this_simname="4s_dave"
this_list=['ba_Ms2.0_Ma0.5_256']
this_list=['aa_Ms2.0_Ma0.5_512']


if 0:
    caq.comp_all(this_list)
if 1:
    p1.plot_quan(this_list)

if 0:
    m2d.make_all(this_list)
if 0:
    maf.make_all(this_list)
if 0:
    m3d.make_spec(this_list)
if 0:
    p3.plot_all_spectra(this_list)

if 0:
    for sim in this_list:
        p1.plot_quan([sim])
if 0:
    imf.image(this_simname)


