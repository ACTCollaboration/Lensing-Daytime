#!/usr/bin/env python
# coding: utf-8

import numpy as np, sys, local, tools_lens, quad_func

kwargs = {'snmin':0,'snmax':100}
kwargs_ov = {'overwrite':True}

kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs)
wn = np.ones(5)
ocl = np.ones((3,aobj_c.lmax+1))
ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
    #for qid in ['02','03','04']:
    for qid in ['diff_dn_234']:
        aobj_d = local.init_analysis_params(qid=qid,ascale=5.0,**kwargs)
        tools_lens.interface(aobj_d,wn,run=['norm','qrec','n0','rdn0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
#tools_lens.interface(aobj_d,wn,run=['norm','qrec','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
#tools_lens.interface(aobj_d,wn,run=['norm','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)


