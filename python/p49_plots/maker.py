from GL import *

import simulation
reload(simulation)
import simulation_info.all_sims as all_sims

import p49_plots.F1_proj as F1
import p49_plots.F2_nominal as F2
import p49_plots.F3_spectra as F3
import p49_plots.F4_amps    as F4
import p49_plots.F5_amp_amp as F5
import p49_plots.F9_bilinear_run as F9
import p49_plots.F6_ratio as F6
import p49_plots.F7_pearson as F7   
import p49_plots.F8_summary as F8
reload(F1);reload(F2);reload(F3)
reload(F4);reload(F5);reload(F6)
reload(F7);reload(F8);reload(F9)

sim_list = all_sims.lists['suite1']


sim_colors.set_cmap('viridis')
for ns, sim in enumerate(sim_list):
    this_sim = simulation.corral[sim]
    this_sim.load()
    #simulation.set_colors(this_sim)
    simulation.set_colors(this_sim,cmap_name=sim_colors.cmap)

proj_cmap = 'plasma'
#proj_cmap = 'cool'
#F1.proj(field='density_',LOS='y', cmap=proj_cmap)
#F1.proj(field='E',LOS='y', no_mean=True, cmap=proj_cmap)
#F1.proj(field='B',LOS='y', cmap=proj_cmap)
#F2.nom(sim_list)
#F3.plot_avg_spectra(sim_list, prim_or_teb='teb', axis='y')
#F3.plot_avg_spectra(sim_list, prim_or_teb='prim')
##F4.plot_amps_slopes(sim_list, prim_or_teb='prim',amps_or_slopes='slopes')
##F4.plot_amps_slopes(sim_list, prim_or_teb='prim',amps_or_slopes='amps')
##F4.plot_amps_slopes(sim_list, prim_or_teb='teb',amps_or_slopes='slopes', axis='y')
##F4.plot_amps_slopes(sim_list, prim_or_teb='teb',amps_or_slopes='amps', axis='y')
#F4.plot_slopes(sim_list, prim_or_teb='prim')
#F4.plot_slopes(sim_list, prim_or_teb='teb')
#F6.plot_ratios(sim_list, LOS='y')
#F6.plot_amps(sim_list,LOS='y')
#F7.plot_spectra(sim_list,LOS='y')
F7.plot_machmean(sim_list,LOS='y')
#F8.plot_summary(sim_list, LOS='y')
#F9.fit_all(sim_list)

if 0:
    fig,ax=plt.subplots(2,1, figsize=(4,8))
    F8.plot_summary(sim_list, LOS='y', ax=ax[1])
    F7.plot_tb_eb(sim_list,LOS='y',ax=ax[0])
    fig.tight_layout()
    fig.savefig("%s/things"%dl.plotdir)

#not part of the main paper
#import F5b_TT_rho as F5b
#reload(F5b)
#F5b.plot_amps_amps(sim_list,amps_or_slopes='slopes',axis='y')

#import F11_gamma as F11
#reload(F11)
#F11.plot_gamma(sim_list)

#import F12_spectratest as F12
#reload(F12)
#F12.plot(['4_1'])#, Nzones2, Nzones1)

#import rework_spectra as rw
#reload(rw)
#rw.rework_spectra(sim_list)


