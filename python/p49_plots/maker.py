from GL import *

import simulation
import simulation_info.all_sims as all_sims

import p49_plots.F2_nominal as F2
import p49_plots.F3_spectra as F3
import p49_plots.F4_amps    as F4
import p49_plots.F5_amp_amp as F5
import p49_plots.F9_bilinear_run as F9
reload(F9)

sim_list = all_sims.lists['suite1']



if 0:
    F2.nom(sim_list)
    F3.plot_avg_spectra(sim_list, prim_or_teb='teb', axis='y')
    F3.plot_avg_spectra(sim_list, prim_or_teb='prim')
    F4.plot_amps_slopes(sim_list, prim_or_teb='prim',amps_or_slopes='slopes')
    F4.plot_amps_slopes(sim_list, prim_or_teb='prim',amps_or_slopes='amps')
    F4.plot_amps_slopes(sim_list, prim_or_teb='teb',amps_or_slopes='slopes', axis='y')
    F4.plot_amps_slopes(sim_list, prim_or_teb='teb',amps_or_slopes='amps', axis='y')
    F5.plot_amps_amps(sim_list, amps_or_slopes='amps',axis='y')
    F5.plot_amps_amps(sim_list, amps_or_slopes='slopes',axis='y')
F9.fit_all(sim_list)
