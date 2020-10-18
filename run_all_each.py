#!/usr/bin/env python
# coding: utf-8

import numpy as np, sys, local, tools_lens, quad_func, tools_cmb

# global fixed parameters
ow = False
vb = False

kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wtype':'iso15v3pt'}

# individual
for qid_d, qid_n in zip(local.boss_d,local.boss_n):
    tools_cmb.diff_day_night(qid_d,qid_n,'diff_'+qid_n,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

kwargs_ov = {'overwrite':ow}
kwargs0 = {'snmin':1,'snmax':100,'ascale':3.0,'wtype':'com16v3'}

kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs0)
ocl = np.ones((3,aobj_c.lmax+1))
ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

for qid_d, qid_n in zip(local.boss_d,local.boss_n):

    for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
        
        aobj_d = local.init_analysis_params(qid='diff_'+qid_n,**kwargs)
        tools_lens.interface(aobj_d,run=['norm','qrec','n0','rdn0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)

