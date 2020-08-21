#!/usr/bin/env python
# coding: utf-8

import numpy as np, curvedsky as cs, sys, plottools as pl, local, tools_cmbalm, tqdm, tools_lens, quad_func
from pixell import enmap

kwargs = {'snmin':0,'snmax':10}
kwargs_ov = {'overwrite':True}
kwargs_qrec = {'n0max':5,'mfmax':10,'qlist':['TT']}

aobj = local.init_analysis_params(qid='diff_dn',**kwargs)
cobj = local.init_analysis_params(qid='comb_d',**kwargs)
wn = np.ones(5)
ocl = np.ones((3,aobj.lmax+1))
ocl[0,:] = (np.loadtxt(aobj.fscl['c'])).T[1]
tools_lens.interface(aobj,wn,run=['norm','qrec','n0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)


