from GL import *

import simulation

import tight_plots
reload(tight_plots)

def proj(field='density_',LOS='y', cmap="winter",no_mean=True):

    is_density=False
    if field=='density_':
        is_density=True

    MACHS = ['half','3','6']
    ALF   = ['half','1','2']


    fig,axes,ccc = tight_plots.fig_squares(len(MACHS),len(ALF))

    array_array=[]
    ext=dt.extents()
    for nm, MMM in enumerate(MACHS):
        for na, AAA in enumerate(ALF):
            name = "%s_%s"%(MMM,AAA)
            this_sim = simulation.corral[name]

            frame = this_sim.ann_frames[-1]
            frb_name = "%s/DD%04d.products/DD%04d_%s%s.fits"%(this_sim.product_location,frame,frame,field,LOS)
            arr = pyfits.open(frb_name)[0].data
            if is_density:
                pass
                #arr=np.log(arr)

            if no_mean:
                arr-=arr.mean()
            array_array.append(arr)
    nplot=0
    def hhh(ttt):
        if ttt == "half":
            Mtext=r"\frac{1}{2}"
            Mtext="0.5"
        else:
            Mtext="%s"%ttt
        return Mtext
    for nm, MMM in enumerate(MACHS):
        for na, AAA in enumerate(ALF):
            name = "%s_%s"%(MMM,AAA)
            label=r'$%s,%s$'%(hhh(MMM),hhh(AAA))

            this_sim = simulation.corral[name]
            norm = mpl.colors.SymLogNorm(linthresh = 0.1, vmin=-1.1,vmax=1.1, base=np.e)
            thax=axes[nm][na]
            arr=array_array[nplot]
            nplot+=1
            plot=thax.imshow(arr, origin='lower',interpolation='nearest',norm=norm, cmap=cmap)
            cmap_tool=mpl.cm.get_cmap(cmap)
            if 0:
                text_color=cmap_tool(256)
                box_color=list(cmap_tool(0))
            else:
                text_color = cmap_tool(0)
                box_color  = list(cmap_tool(256))
            box_color[3]=0.5
            edge_color=cmap_tool(128)
            thax.set(xticks=[],yticks=[])
            thax.text(50,70,label,c=text_color, bbox={'boxstyle':'round','ec':edge_color, 'fc':box_color})
    cb=fig.colorbar(plot, cax=ccc)
    outname='%s/proj_%s.pdf'%(dl.plotdir,field)
    fig.savefig(outname)
    print(outname)




