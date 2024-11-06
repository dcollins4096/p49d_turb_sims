from GL import *

import simulation
import simulation_info.all_sims as all_sims

import plots.P1_plot_quan as p1
import plots.P2_image_all as p2
import plots.P3_all_spectra as p3
import plots.P4_spectra_time as p4
import plots.P5_pdfs as p5
import plots.P6_mean_var as p6
import plots.P7_dt_tool as p7
reload(p1)
reload(p2)
reload(p3)
reload(p4)
reload(p5)
reload(p6)
reload(p7)


sim_list = all_sims.lists['suite3']
#sim_list = ['aa_Ms2.0_Ma0.5_512']
#sim_list = ['4_1']#,'1_1']

if 1:
    p1.plot_all_mach(sim_list)

if 1:
    import simulation as sim
    def mach_arrays(sim_list):
        ms_nom=[]
        ma_nom=[]
        ms_act=[]
        for ns, sim_name in enumerate(sim_list):
            this_sim=sim.corral[sim_name]
            this_sim.read_avg_quan()
            ms_nom.append(this_sim.Ms_nom)
            ma_nom.append(this_sim.Ma_nom)
            QQQ = this_sim.quan_time
            ms_act.append(QQQ['vrms'].mean()/np.sqrt(3))
        fptr = open('%s/mach.txt'%plot_dir,'w')
        for i in range(len(ms_nom)):
            fptr.write('%0.12f %0.12f %0.12f\n'%(ms_nom[i],ma_nom[i],ms_act[i]))
        fptr.close()
        #print("Ms",ms_nom)
        #print("Ma",ma_nom)
        #print("Mr",ms_act)
        return nar(ms_nom),nar(ma_nom),nar(ms_act)
    ms_nom,ma_nom,ms_act=mach_arrays(sim_list)

if 0:
    ok = (ma_nom == 0.0)
    #print(ms_act[ok])
    #interpo = np.interp( ms_nom[ok], ms_act[ok], ms_nom[ok])
    #print(interpo)
    interpo = np.interp( ms_nom[ok], ms_act[ok], ms_nom[ok]/ms_act[ok])*ms_nom[ok]
    print(interpo)





if 0:
    #all plots in one pannel
    p1.plot_quan(sim_list)

if 0:
    #frames can be "all" or "ann"
    p3.plot_all_spectra(sim_list, all_or_ann='ann', compensate=False)

if 0:
    #12 panel image for each frame
    #Plotting all at once puts them all on one 8 
    for sim in sim_list:
        p2.image(sim)

if 0:
    p7.plot_dt(sim_list)

if 0:
    for sim in sim_list:
        this_sim = simulation.corral[sim]
        print(this_sim.ann_frames)

if 0:
    #quan plot, each sim
    for sim in sim_list:
        p1.plot_quan([sim])


if 0:
    p4.slope_time(sim_list)

if 0:
    #PDF of magnetic field components
    sim_list = all_sims.lists['suite1']
    fields=['magnetic_field_%s'%s for s in 'xyz']
    if 0:
        #everything, linear, looks pretty gaussian
        p5.plot_pdfs(sim_list,fields, name = "All_Linear", pdf_prefix='pdf_scaled', all_or_ann_frames='ann', norm_axis=False, overgauss=False, logy=False, plot_all=True)
    if 1:
        #everyting, log, normalized
        p5.plot_pdfs(sim_list,fields, name = "All_log_norm", pdf_prefix='pdf_scaled', all_or_ann_frames='ann', norm_axis=False, overgauss=True, logy=True, plot_all=True)
    if 1:
        #everyting, log, normalized
        p5.plot_pdfs(sim_list,fields, name = "3avg", pdf_prefix='pdf_scaled', all_or_ann_frames='ann', norm_axis=False, overgauss=False, logy=False, plot_all=False, all_on_one=True)
    #p5.plot_pdfs(sim_list,fields, name = "raw_norm", pdf_prefix='pdf', all_or_ann_frames='all',norm_axis=True) #have to use "all" frames with "pdf" prefix

if 0:
    sim_list = all_sims.lists['suite1']
    fields=['magnetic_field_%s'%s for s in 'xyz']
    p5.plot_sigma(sim_list)

if 0:
    #this makes magnetic field PDFs that work.
    #sim_list = ['4_1','1_1','half_half']
    sim_list = all_sims.lists['suite1']
    #sim_list = ['1_1','5_1']
    p5.plot_pdfs_fits(sim_list, name = "Bfields_J", norm_axis=True)

if 0:
    #gets the better answer
    sim_list = all_sims.lists['suite1']
    fields=['magnetic_field_strength']
    p5.plot_pdfs(sim_list,fields, name = "raw_norm", pdf_prefix='pdf', all_or_ann_frames='all',norm_axis=True) #have to use "all" frames with "pdf" prefix

if 0:
    #working here 7/24
    #Use this for all PDFs.
    #sim_list = ['4_1']#,'1_1']
    sim_list = all_sims.lists['suite1']
    fields=['magnetic_field_strength']
    p5.plot_pdfs(sim_list,fields, name = "Bfields", pdf_prefix='pdf_scaled', all_or_ann_frames='ann', norm_axis=False)
    #p5.plot_pdfs(sim_list,fields, name = "raw_norm", pdf_prefix='pdf', all_or_ann_frames='all',norm_axis=True) #have to use "all" frames with "pdf" prefix


if 0:
    #get Btotal to make sense
    fields =['magnetic_field_strength' ]
    p5.plot_pdfs(sim_list,fields, name = "test", pdf_prefix='pdf')
if 0:
    #get Btotal to make sense
    fields =['magnetic_field_x' ]
    p5.plot_pdfs(sim_list,fields, name = "test", pdf_prefix='pdf')

if 0:
    #use this to find sigma_b for each sim
    sim_list = all_sims.lists['suite1']
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z']#,'magnetic_field_strength' ]
    p5.stack_pdfs_raw(sim_list,fields, name = "Bfields_stretched", pdf_prefix='pdf_scaled')
if 0:
    #use this to find sigma_b for each sim
    sim_list = all_sims.lists['suite1']
    fields =['magnetic_field_x','magnetic_field_y','magnetic_field_z']#,'magnetic_field_strength' ]
    p5.pdf_totals(sim_list, name = "tots")

if 0:
    sim_list = all_sims.lists['suite1']
    p6.mean_var(sim_list)
