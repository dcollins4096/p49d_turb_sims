
import get_all_quantities as gaq
reload(gaq)
import sim_colors 
reload(sim_colors)

plotdir = "//home/dccollins/PigPen"
sim_dir = "/data/cb1/Projects/P49_EE_BB/512_frbs/%s"%sim
sim_dir = "/data/cb1/Projects/P49_EE_BB/"
plt.close('all')
if 'quand' not in dir():
    quand = {}
for i, sim in enumerate(sim_colors.simlist):
    if sim in quand:
        continue
    print("=== %s ==="%sim)
    ol = "%s/%s/OutputLog"%(sim_dir,sim)
    quand[sim] = gaq.all_quan_from_outputlog(ol)



fig, axes = plt.subplots(2,3)
ax_list=axes.flatten()
ax_vel = ax_list[0]
ax_valf = ax_list[1]
ax_msms = ax_list[2]
ax_mama = ax_list[3]
ax_mama2 = ax_list[4]
ax_b2v2 = ax_list[5]

vmeans=[]
bmeans=[]
bmeans_also=[]
whatsit=[]
for i, sim in enumerate(sim_colors.simlist):

    quan=quand[sim]
    vx2 = quan['vx_std']**2
    vy2 = quan['vy_std']**2
    vz2 = quan['vz_std']**2
    vrms = np.sqrt(vx2 + vy2 + vz2)

    vax2 = quan['alf_x_std']**2
    vay2 = quan['alf_y_std']**2
    vaz2 = quan['alf_z_std']**2
    varms = np.sqrt(vax2 + vay2 + vaz2)
    print(quan['time'])
    mytime=quan['time']
    mytime -= mytime[0]
    mytime /= mytime[-1]

    ax_vel.plot( mytime, vrms, sim_colors.glyph[sim])
    ax_vel.plot( mytime, [sim_colors.Ms[sim]]*vrms.size, 'k'+sim_colors.linestyle[sim])
    ax_valf.plot( mytime, varms, sim_colors.glyph[sim])
    ax_valf.plot( mytime, [sim_colors.Ma[sim]]*vrms.size, 'k'+sim_colors.linestyle[sim])

    ax_msms.scatter( sim_colors.Ms[sim], vrms.mean(), c=sim_colors.color[sim],marker=sim_colors.marker[sim])
    ax_mama.scatter( sim_colors.Ma[sim], varms.mean(), c=sim_colors.color[sim],marker=sim_colors.marker[sim])

    Bx = quan['bx_avg']
    By = quan['by_avg']
    Bz = quan['bz_avg']
    Btot = np.sqrt( Bx**2+By**2+Bz**2)
    Ma2 = vrms/Btot

    ax_mama2.scatter( sim_colors.Ma[sim], Ma2.mean(), c=sim_colors.color[sim],marker=sim_colors.marker[sim])     

    bx = quan['bx_std']**2
    by = quan['by_std']**2
    bz = quan['bz_std']**2
    brms = np.sqrt( bx + by + bz)
    #ax_b2v2.scatter( vrms.mean(), 1/(vrms.mean()/brms.mean()/Btot.mean()), c=sim_colors.color[sim],marker=sim_colors.marker[sim]) 
    vmeans.append(vrms.mean())
    bmeans.append(brms.mean())
    whatsit.append(Btot.mean()*brms.mean()/vrms.mean())
    bmeans_also.append( vrms.mean()**2/Btot.mean())
    ax_b2v2.scatter( bmeans[-1], bmeans_also[-1], c=sim_colors.color[sim],marker=sim_colors.marker[sim]) 
    #ax_b2v2.scatter( vmeans[-1], whatsit[-1], c=sim_colors.color[sim],marker=sim_colors.marker[sim]) 


pfit=np.polyfit(bmeans, bmeans_also,1)
ax_b2v2.plot( bmeans, 3/2*nar(bmeans), c=[0.5]*3)
ax_msms.plot( sim_colors.ms_list, sim_colors.ms_list, c=[0.5]*4)
ax_mama2.plot( sim_colors.ma_list, sim_colors.ma_list, c=[0.5]*4)
dt.axbonk(ax_vel, xlabel=r'$t$', ylabel=r'$\sigma_{\rm{S}}$',yscale='log')
dt.axbonk(ax_valf, xlabel=r'$t$', ylabel=r'$\sigma_{\rm{A}}$',yscale='log')
dt.axbonk(ax_msms, xlabel=r'$\mathcal{M}_{\rm{S}}$', ylabel=r'$\sigma_{\rm{S}}$')
dt.axbonk(ax_mama, xlabel=r'$\mathcal{M}_{\rm{A}}$', ylabel=r'$\sigma_{\rm{A1}}$')
dt.axbonk(ax_mama2, xlabel=r'$\mathcal{M}_{\rm{A}}$', ylabel=r'$v_{\rm{rms}}/||\langle B\rangle||$')
fig.savefig("%s/quans_time.png"%plotdir)




