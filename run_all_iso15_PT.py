#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb
import tools_lens
import quad_func

# global fixed parameters
kwargs_ov = {'overwrite':True,'verbose':False}
kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wtype':'iso15v3PT','clmin':100}
kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

#run_cmb = ['alm','aps','comb']
run_cmb = ['comb']
#run_cmb = []
#run_cmb = ['alm']

qids_d = ['boss_d03','boss_d04']
qids_n = local.boss_n
qids = ['boss_04']
#qids = qids_d + qids_n

tools_cmb.interface(run_cmb,qids,kwargs_ov=kwargs_ov,kwargs=kwargs)

aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs)
#kwargs0 = {'snmin':1,'snmax':100,'ascale':3.0,'wtype':'com16v3pt'}
#aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs0)
ocl = np.ones((3,aobj_c.lmax+1))
ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

#for qid_d, qid_n in zip(qids_d,qids_n):

for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
    
    aobj_d = local.init_analysis_params(qid='diff_dn',**kwargs)
    #aobj_d = local.init_analysis_params(qid='diff_'+qid_n,**kwargs)
    tools_lens.interface(aobj_d,run=['norm','qrec','n0','rdn0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
    #tools_lens.interface(aobj_d,run=['norm','qrec'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)

