#!/usr/bin/env python
# coding: utf-8

import numpy as np, sys, os
import local
import tools_lens
import quad_func

kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wind':'com16','ivar':'v0'}
kwargs_ov = {'overwrite':False}

kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
#for kwargs_qrec in [kwargs_qrec0]:
    #for qid in ['02','03','04']:
    for qid in ['comb_d','comb_n']:

        aobj_c = local.init_analysis_params(qid=qid,**kwargs)
        ocl = np.ones((3,aobj_c.lmax+1))
        ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

        aobj_d = local.init_analysis_params(qid=qid,**kwargs)
        #tools_lens.interface(aobj_d,wn,run=['rdn0'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
        tools_lens.interface(aobj_d,run=['norm','qrec','n0','rdn0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
#tools_lens.interface(aobj_d,wn,run=['norm','qrec','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
#tools_lens.interface(aobj_d,wn,run=['norm','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)

