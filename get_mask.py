#!/usr/bin/env python
# coding: utf-8

import numpy as np, healpy as hp, curvedsky as cs, sys
import local, tools_cmb

sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/actsims/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/soapack/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/orphics/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/tilec/")
from pixell import enmap
from soapack import interfaces


qids = ['boss_01','boss_d02','boss_02','boss_d03','boss_03','boss_d04','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['boss_d01','boss_01','boss_d02','boss_02','boss_d03','boss_03','boss_d04','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['boss_d01']

for ap in [1.,2.,3.,4.,5.]:

    aobj = {q: local.init_analysis_params(qid=q,ascale=ap) for q in qids}

    for q in qids:
        
        mask_iv = tools_cmb.load_mask(q,with_ivar=False)
        mask = enmap.to_healpix(mask_iv[0],nside=2048)
    
        if ap != 1.:
            mask_binary = mask/(mask+1e-30)
            print(q,np.average(mask_binary))
            amask = cs.utils.apodize(2048,mask_binary,aobj[q].ascale)
        else:
            amask = mask.copy()
        
        print(q,np.average(amask))        
        hp.fitsfunc.write_map(aobj[q].amask,amask,overwrite=True)
