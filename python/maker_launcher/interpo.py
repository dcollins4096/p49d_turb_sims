#!/usr/bin/env python
import numpy as np

fptr = open('mach.txt','r')
lines = fptr.readlines()
fptr.close()
ms_nom=[]
ma_nom=[]
ms_act=[]
for line in lines:
    tok = line[:-1].split()
    ms_nom.append(float(tok[0]))
    ma_nom.append(float(tok[1]))
    ms_act.append(float(tok[2]))

ms_nom=np.array(ms_nom)
ma_nom=np.array(ma_nom)
ms_act=np.array(ms_act)

def get_ms_nom(ms,ma):
    ok = ma_nom == ma
    if ok.sum() == 0 :
        print("Invalid ma", ma, "not in ", np.unique(ma_nom))
    else:
        msn = ms_nom[ok]
        msa = ms_act[ok]
        ms_choice = ms*np.interp(ms, msa, msn/msa)
        return(ms_choice)

#end
