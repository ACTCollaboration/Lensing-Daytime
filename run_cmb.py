#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb

# global fixed parameters
ow = True
vb = False

qids = local.qid_all
#kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wtype':'com16'}
kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wtype':'com15'}
#kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wtype':'base'}

#run_cmb = ['alm','aps','comb']
#run_cmb = ['map','alm','aps','sup']
#run_cmb = ['alm','aps','sup']
#run_cmb = ['alm','aps','each']
run_cmb = ['each']


if 'map' in run_cmb:
    print('Generating map')
    tools_cmb.generate_map(qids,overwrite=ow,verbose=vb,**kwargs)

if 'alm' in run_cmb:
    print('Convert map to alm')
    tools_cmb.map2alm(qids,overwrite=ow,**kwargs)

if 'aps' in run_cmb:
    print('Compute aps')
    tools_cmb.alm2aps(qids,overwrite=ow,verbose=vb,mtype=['T'],**kwargs)

if 'sup' in run_cmb:
    print('Supfuc')
    tools_cmb.alm_supfac(qids,overwrite=ow,verbose=vb,**kwargs)


if 'each' in run_cmb:

    # individual
    for qid_d, qid_n in zip(local.boss_d,local.boss_n):
        qid_diff = 'diff_'+qid_n
        # day - night
        tools_cmb.diff_day_night(qid_d,qid_n,qid_diff,overwrite=ow,verbose=True,mtype=['T'],**kwargs)


if 'comb' in run_cmb:
    
    #qids_d = local.boss_d
    #qids_n = local.boss_n
    #qid_dn = ['comb_d','comb_n']
    #qid_diff = 'diff_dn'

    #qids_d = local.boss_d[1:]
    #qids_n = local.boss_n[1:]
    #qid_dn = ['comb_d_234','comb_n_234']
    #qid_diff = 'diff_dn_234'

    if wtype == 'base':
        sys.exit('not support for wtype=base')
    
    if '15' in wtype:
        qids_d = local.boss_d
    else:
        qids_d = local.boss_d + local.s_16_d

    qids_n = local.boss_n
    qid_dn = ['comb_d','comb_n']
    qid_diff = 'diff_dn'

    # comb day or night
    print('combine alm for day')
    tools_cmb.alm_comb(qids_d,qid_dn[0],overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    print('combine alm for night')
    tools_cmb.alm_comb(qids_n,qid_dn[1],overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    
    # W2 factors from each day or night masks
    print('compute aps')
    tools_cmb.alm2aps(qid_dn,overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    #print('compute supfac')
    #tools_cmb.alm_supfac(qid_dn,W1={'comb_d':0.016,'comb_n':0.02},overwrite=ow,verbose=False,**kwargs)

    # day - night
    print('compute day - night')
    tools_cmb.diff_day_night(qid_dn[0],qid_dn[1],qid_diff,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

    # day + night
    tools_cmb.alm_comb(qids_d+qids_n,'comb_dn',overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    tools_cmb.alm2aps(['comb_dn'],overwrite=ow,verbose=False,mtype=['T'],W2=0.009,**kwargs)

