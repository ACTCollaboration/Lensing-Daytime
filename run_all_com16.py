#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb
import tools_lens
import quad_func

# global fixed parameters
kwargs_ov = {'overwrite':False,'verbose':False}
kwargs = {'snmin':30,'snmax':100,'ascale':3.0,'wtype':'com16v0pt','clmin':100}
kwargs_qrec0 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':50,'mfmax':100,'rlmin':500,'qlist':['TT'],'rd4sim':True}

#run_cmb = ['comb']
#run_cmb = ['alm','aps','comb']
run_cmb = []
qids = local.qid_all

tools_cmb.interface(run_cmb,qids,kwargs_ov=kwargs_ov,kwargs=kwargs)

#for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
for kwargs_qrec in [kwargs_qrec1]:
    #for qid in ['comb_dn','comb_n']:
    for qid in ['comb_n']:
        aobj = local.init_analysis_params(qid=qid,**kwargs)
        #tools_lens.interface(aobj,run=['norm','qrec','n0','rdn0','mean','aps'],kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
        tools_lens.interface(aobj,run=['rdn0'],kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)

