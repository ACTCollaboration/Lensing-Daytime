#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb
import tools_lens
import quad_func

# global fixed parameters
kwargs_ov = {'overwrite':True,'verbose':False}
kwargs = {'snmin':0,'snmax':100,'ascale':3.0,'wind':'com15','ivar':'noivar','clmin':100}
kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT']}

#run_cmb = ['comb']
run_cmb = ['alm','aps','diff']
#run_cmb = []
qids = local.boss_dn

tools_cmb.interface(run_cmb,qids,kwargs_ov=kwargs_ov,kwargs=kwargs)

for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
    for qid in ['comb_dn','comb_d']:
    #for qid in ['comb_n']:
        tools_lens.interface(qid,run=['norm','qrec','n0','rdn0','mean','aps'],kwargs_ov=kwargs_ov,kwargs_cmb=kwargs,kwargs_qrec=kwargs_qrec)

