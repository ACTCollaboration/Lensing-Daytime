#!/usr/bin/env python
# coding: utf-8

import numpy as np, sys, local, tools_cmb, tools_lens, quad_func

#kwargs = {'snmin':0,'snmax':1}
kwargs = {'snmin':0,'snmax':100}
kwargs_ov = {'overwrite':True}

kwargs_qrec0 = {'n0max':int(kwargs['snmax']/2),'mfmax':kwargs['snmax'],'rlmin':500,'qlist':['TT'],'bhe':['src']}
kwargs_qrec1 = {'n0max':int(kwargs['snmax']/2),'mfmax':kwargs['snmax'],'rlmin':500,'qlist':['TT']}

aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs)
wn = np.ones(5)
ocl = np.ones((3,aobj_c.lmax+1))
ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

#for qid in ['02','04']:
for qid in ['01','02','03','04']:

    #tools_cmb.diff_day_night('boss_d'+qid,'boss_'+qid,'diff_boss_'+qid,overwrite=True,verbose=True,mtype=['T'],**kwargs)
    #qid_out = 'diff_boss_'+qid
    
    # for overlap with s16
    qid_out = 'diff_boss_'+qid+'_ov'
    tools_cmb.diff_day_night('boss_d'+qid,'boss_'+qid,qid_out,overwrite=True,verbose=True,mtype=['T'],ascale=5.,**kwargs)

    for kwargs_qrec in [kwargs_qrec0,kwargs_qrec1]:
        
        aobj_d = local.init_analysis_params(qid=qid_out,ascale=5.0,**kwargs)
        
        #tools_lens.interface(aobj_d,wn,run=['norm','qrec','aps'],meansub=False,ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
        tools_lens.interface(aobj_d,wn,run=['norm','qrec','n0','mean','aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
        tools_lens.interface(aobj_d,wn,run=['rdn0'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)
        #tools_lens.interface(aobj_d,wn,run=['aps'],ocl=ocl,kwargs_ov=kwargs_ov,kwargs_qrec=kwargs_qrec)



