#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb
import tools_lens
import quad_func

# global fixed parameters
kwargs_ov = {'overwrite':False,'verbose':False}
kwargs = {'snmin':0,'snmax':100,'ascale':1.0,'wtype':'com16v0PT','clmin':100}
kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

#run_cmb = ['comb']
run_cmb = ['alm','aps','comb']
#run_cmb = []
qids = local.qid_all

tools_cmb.interface(run_cmb,qids,kwargs_ov=kwargs_ov,kwargs=kwargs)

#for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
for kwargs_qrec in [kwargs_qrec1]:
    for qid in ['comb_dn','comb_d','comb_n']:
        tools_lens.interface(qid,run=['norm','qrec','n0','rdn0','mean','aps'],kwargs_ov=kwargs_ov,kwargs_cmb=kwargs,kwargs_qrec=kwargs_qrec)

