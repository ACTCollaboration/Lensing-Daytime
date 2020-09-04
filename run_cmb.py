#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import tools_cmb


# In[2]:

ow = True

qids = ['boss_d01','boss_d02','boss_d03','boss_d04','boss_01','boss_02','boss_03','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['boss_01','boss_02','boss_03','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['s16_d02','s16_d03']
#qids = ['boss_03','boss_04']
#qids = ['bndn_01','bndn_02']
#qids = ['boss_d01']

kwargs = {'snmin':0,'snmax':100,'ascale':5.0}

#run_cmb = ['alm','aps','sup']
#run_cmb = ['map','alm','aps','sup']
#run_cmb = ['each']
run_cmb = ['comb']


if 'map' in run_cmb:
    print('Generating map')
    tools_cmb.generate_map(qids,overwrite=ow,verbose=True,**kwargs)

if 'alm' in run_cmb:
    print('Convert map to alm')
    tools_cmb.map2alm(qids,overwrite=ow,**kwargs)

if 'aps' in run_cmb:
    print('Compute aps')
    tools_cmb.alm2aps(qids,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

if 'sup' in run_cmb:
    print('Supfuc')
    tools_cmb.alm_supfac(qids,overwrite=ow,verbose=True,**kwargs)


if 'each' in run_cmb:

    # individual
    qids_d = ['boss_d01','boss_d02','boss_d03','boss_d04']
    qids_n = ['boss_01','boss_02','boss_03','boss_04']
    
    for qid_d, qid_n in zip(qids_d,qids_n):
        #qid_diff = 'diff_boss_03_raw'
        qid_diff = 'diff_'+qid_n
        # day - night
        tools_cmb.diff_day_night(qid_d,qid_n,qid_diff,overwrite=ow,verbose=True,mtype=['T'],**kwargs)


if 'comb' in run_cmb:
    
    qids_d = ['boss_d01','boss_d02','boss_d03','boss_d04']
    qids_n = ['boss_01','boss_02','boss_03','boss_04']
    qid_dn = ['comb_d','comb_n']
    qid_diff = 'diff_dn'

    #qids_d = ['boss_d02','boss_d03','boss_d04']
    #qids_n = ['boss_02','boss_03','boss_04']
    #qid_dn = ['comb_d_234','comb_n_234']
    #qid_diff = 'diff_dn_234'
    
    # comb day or night
    print('combine alm for day')
    tools_cmb.alm_comb(qids_d,qid_dn[0],overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    print('combine alm for night')
    tools_cmb.alm_comb(qids_n,qid_dn[1],overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    
    # W2 factors from each day or night masks
    print('compute aps')
    tools_cmb.alm2aps(qid_dn,overwrite=ow,verbose=False,mtype=['T'],W2={'comb_d':0.007,'comb_n':0.01},**kwargs)
    print('compute supfac')
    tools_cmb.alm_supfac(qid_dn,W1={'comb_d':0.016,'comb_n':0.02},overwrite=ow,verbose=False,**kwargs)

    # day - night
    print('compute day - night')
    tools_cmb.diff_day_night(qid_dn[0],qid_dn[1],qid_diff,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

    # day + night
    #tools_cmb.alm_comb(qids_d+qids_n,'comb_dn',overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    #tools_cmb.alm2aps(['comb_dn'],overwrite=ow,verbose=False,mtype=['T'],W2=0.009,**kwargs)

