#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb
import tools_lens
import quad_func

# global fixed parameters
kwargs_ov = {'overwrite':False,'verbose':False}
kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wtype':'iso15v3pt'}
kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

#run_cmb = ['alm','aps','sup','each','comb']
run_cmb = ['alm','aps','comb']

tools_cmb.interface(run_cmb,local.boss_dn,kwargs_ov=kwargs_ov,kwargs=kwargs)

aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs)
ocl = np.ones((3,aobj_c.lmax+1))
ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
    for q in ['boss_02','boss_03','boss_04']:
        aobj_d = local.init_analysis_params(qid='diff_'+q,**kwargs)
        tools_lens.interface(aobj_d,run=['norm','qrec','n0','rdn0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)



'''
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
    for qid_d, qid_n in zip(local.boss_d,local.boss_n):
        tools_cmb.diff_day_night(qid_d,qid_n,'diff_'+qid_n,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

if 'comb' in run_cmb:
    
    qid_dn = ['comb_d','comb_n']

    # comb day or night
    print('combine alm for day')
    tools_cmb.alm_comb(local.boss_d,qid_dn[0],overwrite=ow,verbose=vb,mtype=['T'],**kwargs)
    print('combine alm for night')
    tools_cmb.alm_comb(local.boss_n,qid_dn[1],overwrite=ow,verbose=vb,mtype=['T'],**kwargs)
    
    # W2 factors from each day or night masks
    print('compute aps')
    tools_cmb.alm2aps(qid_dn,overwrite=ow,verbose=vb,mtype=['T'],**kwargs)
    print('compute supfac')
    tools_cmb.alm_supfac(qid_dn,overwrite=ow,verbose=vb,**kwargs)

    # day - night
    print('compute day - night')
    tools_cmb.diff_day_night(qid_dn[0],qid_dn[1],'diff_dn',overwrite=ow,verbose=vb,mtype=['T'],**kwargs)

    # day + night
    tools_cmb.alm_comb(local.boss_dn,'comb_dn',overwrite=ow,verbose=vb,mtype=['T'],**kwargs)
    tools_cmb.alm2aps(['comb_dn'],overwrite=ow,verbose=vb,mtype=['T'],**kwargs)

'''

