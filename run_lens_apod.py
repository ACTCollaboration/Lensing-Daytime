#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_lens
import quad_func

kwargs_cmb = {\
    'snmin':0, \
    'snmax':100, \
    'fltr':'none', \
    'ivar':'noivar' \
}

kwargs_ov = {\
    'overwrite':True, \
    'verbose':True \
}

kwargs_q0 = {\
    'rlmin':500, \
    'qlist':['TT'], \
    'n0max':int(kwargs_cmb['snmax']/2), \
    'mfmax':kwargs_cmb['snmax'], \
    'rdmax':kwargs_cmb['snmax'] \
}

kwargs_q1 = kwargs_q0.copy()
kwargs_q1['bhe'] = ['src']

for kwargs_q in [kwargs_q0,kwargs_q1]:
#for kwargs_q in [kwargs_q0]:
    #for qid in ['boss_s15n']:
    #for qid in ['boss_s15dn']:
    #for qid in ['boss_s16d']:
    for qid in ['boss_alldn']:
        tools_lens.interface(qid,run=['norm','qrec','n0','rdn0','mean','aps'],kwargs_ov=kwargs_ov,kwargs_cmb=kwargs_cmb,kwargs_qrec=kwargs_q)

