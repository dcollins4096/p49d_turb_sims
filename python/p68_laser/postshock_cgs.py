from dtools.starter1 import *
import dtools.davetools as dt
import preshock_regions as preshock
#reload(preshock)
import physical_values as vals
#reload(vals)
plt.close('all')
list_of_shots=['r60_t2']
for shot in list_of_shots:

    I_pre=preshock.denoise_preshock[shot]
    I0_code = I_pre.max()

    I_0 = vals.compute_I0(I0_code)
    print(I_0)
    rho_pre = vals.compute_rho(preshock.preshock_region[shot], I_0)
    rho_pre=rho_pre.real
    rho_post = vals.compute_rho(preshock.postshock_region[shot], I_0)
    rho_post=rho_post.real


    fig,axes=plt.subplots(2,2,figsize=(12,12))
    ax0=axes[0][0]; ax1=axes[0][1];ax2=axes[1][0];ax3=axes[1][1]
    ext=dt.extents()
    ext(rho_pre)
    norm_pre = mpl.colors.Normalize(ext.minmax[0],ext.minmax[1])
    ext(rho_post)
    norm_post = mpl.colors.Normalize(ext.minmax[0],ext.minmax[1])
    plot=ax0.imshow(rho_pre,norm=norm_pre)
    fig.colorbar(plot,ax=ax0)
    plot=ax1.imshow(rho_post,norm=norm_post)
    fig.colorbar(plot,ax=ax1)





    fig.savefig('plots_to_sort/%s_cgs'%shot)







