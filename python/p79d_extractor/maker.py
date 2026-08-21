
from  GL import *

import simulation
import simulation_info.all_sims as all_sims
import simulation_info.suite_2

import p79d_extractor.puller as puller
reload(puller)
half=None

if 0:
    sim_list = all_sims.lists['suite1b']
    size=512
    N_per_frame = 5
    rotate = False
    target_res = 128
    suffix='_suite1b_test'
if 0:
    sim_list = []
    for sim in all_sims.lists['suite1']:
        if sim[0] in '2356':
            sim_list.append(sim)
    size=512
    N_per_frame = 5
    rotate = False
    target_res = 128
    los = 'xyz'
    suffix='_2356'
    print(sim_list)
if 0:
    sim_list = []
    for sim in all_sims.lists['suite1']:
        if sim[0] in '4':
            sim_list.append(sim)
    size=512
    N_per_frame = 5
    rotate = False
    target_res = 128
    los = 'xyz'
    suffix='_4'
    print(sim_list)
if 0:
    sim_list = []
    for sim in all_sims.lists['suite1']:
        if sim[0] in '23456':
            sim_list.append(sim)
    size=512
    N_per_frame=5
    rotate=False
    target_res=128
    los='xyz'
    half=1
    suffix='23456'
if 0:
    sim_list = ['5_1']
    size=512
    N_per_frame=5
    rotate=False
    target_res=128
    los='y'
    half=0
    suffix='_5-1'
if 0:
    sim_list = []
    for sim in all_sims.lists['suite1']:
        if sim[0] in '23456':
            sim_list.append(sim)
    size=256
    N_per_frame=10
    rotate=True
    target_res=128
    los='xyz'
    half=1
    suffix='23456'
if 0:
    sim_list = []
    for sim in all_sims.lists['suite1']:
        if sim[0] in '23456':
            sim_list.append(sim)
    size=256
    N_per_frame=5
    rotate=False
    target_res=128
    los='xyz'
    half=1
    suffix='23456'

if 0:
    sim_list = []
    for sim in all_sims.lists['suite1']:
        if sim[0] in '23456':
            sim_list.append(sim)
    size=256
    N_per_frame=5
    rotate=False
    target_res=None
    los='xyz'
    half=0
    suffix='23456'
if 0:
    #sim_list = all_sims.lists['suite4']
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite4']):
        if ns%4!=3:
            sim_list.append(sim)
    size=256
    N_per_frame=5
    rotate=False
    target_res=64
    los='xyz'
    half=1
    suffix='suite4_QU_'
if 0:
    #sim_list = all_sims.lists['suite4']
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite4']):
        if ns%4!=3:
            sim_list.append(sim)
    size=256
    N_per_frame=5
    rotate=False
    target_res=64
    los='y'
    half=0
    suffix='fixed_mach_QU'
if 0:
    #sim_list = all_sims.lists['suite4']
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite4']):
        if ns%4!=3:
            sim_list.append(sim)
    size=512
    N_per_frame=5
    rotate=False
    target_res=64
    los='y'
    half=1
    suffix='THQUEB'

fields='THQUEB'
if 0:
    #sim_list = all_sims.lists['suite4']
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite4']):
        if ns%4==3:
            sim_list.append(sim)
    size=512
    N_per_frame=7
    rotate=False
    target_res=128
    los='xyz'
    half=1
    suffix='T_annfix'
    fields='T'
if 0:
    #sim_list = all_sims.lists['suite4']
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite4']):
        if ns%4!=3:
            sim_list.append(sim)
    size=512
    N_per_frame=7
    rotate=False
    target_res=128
    los='xyz'
    half=0
    suffix='T'
    fields='T'
if 0:
    sim_list = all_sims.lists['suite5']
    size=256
    N_per_frame=1
    rotate=False
    target_res=64
    los='xyz'
    half=1
    suffix='suite5'
    fields='T'
if 0:
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite5']):
        if ns <= 9:
            sim_list.append(sim)
    size=256
    N_per_frame=1
    rotate=False
    target_res=64
    los='xyz'
    half=1
    suffix='suite5_machLE5'
    fields='T'
if 0:
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite7']):
        sim_list.append(sim)
    size=128
    N_per_frame=5
    rotate=False
    target_res=None
    los='xyz'
    half=0
    suffix='suite7c'
    fields='T'
if 0:
    sim_list=[]
    for ns,sim in enumerate(all_sims.lists['suite7']):
        sim_list.append(sim)
    size=128
    N_per_frame=1
    rotate=False
    target_res=None
    los='xyz'
    half=0
    suffix='suite7vs'
    fields='TVS'
if 0:
    sim_list=all_sims.lists['suite1']
    #for ns,sim in enumerate(all_sims.lists['suite1']):
    #    if sim[0] not in ['4','5','6']:
    #        continue
    #    sim_list.append(sim)
    size=128
    N_per_frame=5
    rotate=False
    target_res=128
    los='xyz'
    half=1
    suffix='suite1_tvsquhp'
    fields='TVSQUHP'
if 0:
    sim_list=all_sims.lists['brano']
    #for ns,sim in enumerate(all_sims.lists['suite1']):
    #    if sim[0] not in ['4','5','6']:
    #        continue
    #    sim_list.append(sim)
    size=512
    N_per_frame=5
    rotate=False
    target_res=128
    los='xyz'
    half=1
    suffix='brano_tvsquhp'
    fields='TVSQUHP'
if 1:
    sim_list=all_sims.lists['mach_grid'][1:]
    #for ns,sim in enumerate(all_sims.lists['suite1']):
    #    if sim[0] not in ['4','5','6']:
    #        continue
    #    sim_list.append(sim)
    size=256
    N_per_frame=5
    rotate=False
    target_res=128
    los='xyz'
    half=1
    suffix='mach_grid_tvs_missingOne'
    fields='TVS'
print(sim_list)
if 1:
    puller.pull(sim_list, size, N_per_frame, target_res = target_res, suffix=suffix, rotate=rotate, los=los,half=0, fields=fields)
    puller.pull(sim_list, size, N_per_frame, target_res = target_res, suffix=suffix, rotate=rotate, los=los,half=1, fields=fields)
