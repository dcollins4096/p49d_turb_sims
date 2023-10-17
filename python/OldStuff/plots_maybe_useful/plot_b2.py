
import get_all_quantities as gaq
reload(gaq)
import sim_colors 
reload(sim_colors)

plotdir = "//home/dccollins/PigPen"
sim_dir = "/data/cb1/Projects/P49_EE_BB/512_frbs/%s"%sim
sim_dir = "/data/cb1/Projects/P49_EE_BB/"
plt.close('all')


fig, axes = plt.subplots(1,1)
axb2=axes
#ax_list=axes.flatten()
#ax_vel = ax_list[0]
#ax_valf = ax_list[1]
#ax_msms = ax_list[2]
#ax_mama = ax_list[3]
#ax_mama2 = ax_list[4]
#ax_b2v2 = ax_list[5]

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




    Bx = quan['bx_avg']
    By = quan['by_avg']
    Bz = quan['bz_avg']
    Btot = np.sqrt( Bx**2+By**2+Bz**2)
    Ma2 = vrms/Btot


    bx = quan['bx_std']**2
    by = quan['by_std']**2
    bz = quan['bz_std']**2
    brms = np.sqrt( bx + by + bz)
    #ax_b2v2.scatter( vrms.mean(), 1/(vrms.mean()/brms.mean()/Btot.mean()), c=sim_colors.color[sim],marker=sim_colors.marker[sim]) 
    vmeans.append(vrms.mean())
    bmeans.append(brms.mean())
    whatsit.append(Btot.mean()*brms.mean()/vrms.mean())
    bmeans_also.append( vrms.mean()**2/Btot.mean())

    #axb2.plot( mytime, 1.5*brms, sim_colors.glyph[sim])
    axb2.plot( mytime, vrms**2/Btot/(1.5*brms), sim_colors.glyph[sim])


dt.axbonk(axb2, xlabel=r'$t$', ylabel=r'$v_{\rm{rms}}/||\langle B\rangle||$')
fig.savefig("%s/b2_time.png"%plotdir)




