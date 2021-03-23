#!/usr/bin/env python
# coding: utf-8

import numpy as np, sys, os
import local
import tools_lens
import quad_func

kwargs = {'snmin':1,'snmax':50,'ascale':3.0,'wtype':'com16v3pt'}
kwargs_ov = {'overwrite':False}
kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'rd4sim':True}

for kwargs_qrec in [kwargs_qrec0]:
    for qid in ['comb_dn']:
        aobj_d = local.init_analysis_params(qid=qid,**kwargs)
        tools_lens.interface(aobj_d,run=['rdn0'],kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
