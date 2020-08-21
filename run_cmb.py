#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import tools_cmbalm


# In[2]:

ow = True

qids = ['boss_d01','boss_d02','boss_d03','boss_d04','boss_01','boss_02','boss_03','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['boss_01','boss_02','boss_03','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['s16_d02','s16_d03']
#qids = ['boss_03','boss_04']
#qids = ['bndn_01','bndn_02']
#qids = ['boss_01']

kwargs = {'snmin':0,'snmax':3}

#run_cmb = ['map','alm','aps','sup','dif']
#run_cmb = ['map','alm','aps','sup']
run_cmb = ['comb']


if 'map' in run_cmb:
    print('Generating map')
    tools_cmbalm.generate_map(qids,overwrite=ow,verbose=True,**kwargs)

if 'alm' in run_cmb:
    print('Convert map to alm')
    tools_cmbalm.map2alm(qids,overwrite=ow,**kwargs)

if 'aps' in run_cmb:
    print('Compute aps')
    tools_cmbalm.alm2aps(qids,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

if 'sup' in run_cmb:
    print('Supfuc')
    tools_cmbalm.alm_supfac(qids,overwrite=ow,verbose=True,**kwargs)

if 'comb' in run_cmb:
    qids_d = ['boss_d01','boss_d02','boss_d03','boss_d04']
    qids_n = ['boss_01','boss_02','boss_03','boss_04']
    qid_dn = ['comb_d','comb_n']
    
    # comb day or night
    tools_cmbalm.alm_comb(qids_d,qid_dn[0],overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    tools_cmbalm.alm_comb(qids_n,qid_dn[1],overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    tools_cmbalm.alm2aps(qid_dn,overwrite=ow,verbose=False,mtype=['T'],Wn=np.ones(5),**kwargs)
    tools_cmbalm.alm_supfac(qid_dn,w1=1.,overwrite=ow,verbose=False,**kwargs)
    
    # day - night
    #tools_cmbalm.diff_day_night(overwrite=ow,verbose=True,mtype=['T'],**kwargs)
    
    # day + night
    #tools_cmbalm.alm_comb(qids_d+qids_n,'comb_dn',overwrite=ow,verbose=False,mtype=['T'],**kwargs)
    #tools_cmbalm.alm2aps(['comb_dn'],overwrite=ow,verbose=False,mtype=['T'],Wn=np.ones(5),**kwargs)


