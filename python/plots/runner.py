from GL import *

import simulation
import simulation_info.all_sims as all_sims

import plots.P1_plot_quan as p1
import plots.P2_image_all as p2
import plots.P3_all_spectra as p3
import plots.P4_spectra_time as p4
import plots.P5_pdfs as p5
import plots.P6_mean_var as p6
reload(p1)
reload(p2)
reload(p3)
reload(p4)
reload(p5)
reload(p6)

#sim_list = all_sims.lists['suite1b']
#sim_list = ['4_1','1_1']


if 0:
    for sim in sim_list:
        this_sim = simulation.corral[sim]
        print(this_sim.ann_frames)

if 0:
    #all plots in one pannel
    p1.plot_quan(sim_list)
if 0:
    #quan plot, each sim
    for sim in sim_list:
        p1.plot_quan([sim])
if 0:
    #12 panel image for each frame
    for sim in sim_list:
        p2.image(sim)

if 0:
    #frames can be "all" or "ann"
    p3.plot_all_spectra(sim_list, all_or_ann='ann')

if 0:
    p4.slope_time(sim_list)


if 0:
    #working here bozo
    #sim_list = ['4_1']#,'1_1']
    sim_list = all_sims.lists['suite1']
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength' ]
    #p5.plot_pdfs(sim_list,fields, name = "Bfields", pdf_prefix='pdf_test')
    #p5.plot_pdfs(sim_list,fields, name = "Bfields", pdf_prefix='pdf_scaled')

if 1:
    #bozo here first
    fields =['magnetic_field_strength' ]
    p5.plot_pdfs(sim_list,fields, name = "test", pdf_prefix='pdf')
if 0:
    #use this to find sigma_b for each sim
    sim_list = all_sims.lists['suite1']
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z']#,'magnetic_field_strength' ]
    p5.stack_pdfs_raw(sim_list,fields, name = "Bfields_stretched")
if 0:
    #use this to find sigma_b for each sim
    sim_list = all_sims.lists['suite1']
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z']#,'magnetic_field_strength' ]
    p5.pdf_totals(sim_list, name = "tots")

if 0:
    #this makes magnetic field PDFs that work.
    #sim_list = ['4_1','1_1','half_half']
    sim_list = all_sims.lists['suite1']
    #sim_list = ['1_1','5_1']
    
    p5.plot_pdfs_fits(sim_list, name = "Bfields_J")

if 0:
    sim_list = all_sims.lists['suite1']
    p6.mean_var(sim_list)
