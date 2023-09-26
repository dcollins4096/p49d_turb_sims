from GL import *
import simulation

import itertools

import bilinear
reload(bilinear)


#save bi-linear fit parameters to an h5 file
#function recording values for yaxis (perp to b-field)
def fit_all(simlist,do_plot=False, fit34=3):
    herd = {}
    #field_list = ['avg_clee','avg_clbb','avg_cltt','avg_v','avg_d','avg_h']
    field_list =  ['density','velocity','magnetic','ClTTy','ClEEy','ClBBy']
    #field_list = ['avg_d']
    #field_list = ['avg_clee','avg_clbb','avg_cltt']
    for nf,field in enumerate(field_list):
        msarr=[]
        maarr=[]
        amparr=[]
        slopearr=[]

        for sim in simlist:
            this_sim = simulation.corral[sim]
            this_sim.load()
            if 0:
                msarr=np.append(msarr,this_sim.Ms_mean)
                maarr=np.append(maarr,this_sim.Ma_mean)
            if 1:
                msarr=np.append(msarr,this_sim.quan_mean['msavg'])
                maarr=np.append(maarr,this_sim.quan_mean['maavg'])
            amparr=np.append(amparr,this_sim.ampsA[field])
            slopearr=np.append(slopearr,this_sim.slopesA[field])
        herd[field + "_s"] = bilinear.beefitter(field, msarr, maarr, slopearr, fit34=fit34)
        herd[field + "_a"] = bilinear.beefitter(field, msarr, maarr, np.log(amparr),fit34=fit34)

    field_order =  ['density_s','velocity_s','magnetic_s','ClTTy_s','ClEEy_s','ClBBy_s',
                    'density_a','velocity_a','magnetic_a','ClTTy_a','ClEEy_a','ClBBy_a']
    bilinear.write_tex( herd, '%s/table1.tex'%dl.plotdir, field_order)

    #predict
    truth = {'ClEEy_s':-2.4, 'ClBBy_s':-2.5, 'ClTTy_s':-2.6}
    truth = {'ClEEy_s':-2.4, 'ClBBy_s':-2.5}#, 'ClTTy_s':-2.6}
    truth_fields=truth.keys()

    print('real values')
    pred_ms = []; pred_ma = []
    for A, B in itertools.combinations(truth_fields,2):
        predict_Ms, predict_Ma = bilinear.pairwise(herd[A], herd[B], truth[A],truth[B])
        print("Predict: Ms %0.3f Ma %0.3f"%(predict_Ms,predict_Ma))
        pred_ms.append(predict_Ms)
        pred_ma.append(predict_Ma)
        reload(bilinear)
        #LLL=[herd[A]]
        if do_plot:
            LLL=[herd[A], herd[B]]
            bilinear.plot_herd(LLL, truth=[truth[A],truth[B]], predict=[predict_Ms,predict_Ma])
    return herd


