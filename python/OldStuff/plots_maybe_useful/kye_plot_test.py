
fig2,ax2 = plt.subplots(1,1)
for i in range(len(simdic)):
    ax2.scatter( alfmachavg[i],sonicmachavg[i])
fig2.savefig('%s/msmameas.png'%plotdir)
