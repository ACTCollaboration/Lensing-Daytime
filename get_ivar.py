#!/usr/bin/env python
# coding: utf-8

import numpy as np, healpy as hp, curvedsky as cs, sys, local, tools_cmb
from pixell import enmap


for wtype in ['comdy']:
#for wtype in ['com15','com16']:
    
    aobj = local.init_analysis_params(qid='boss_d03',ascale=3.,wtype=wtype)

    if wtype in ['com15','com16']:
        amask = hp.fitsfunc.read_map(aobj.amask)
        amask[amask!=0] = 1.
    else:
        amask = 1.

    if wtype == 'com16':
        qids = local.qid_all
        fivar = aobj.fivar16
    if wtype == 'com15':
        qids = local.boss_dn
        fivar = aobj.fivar15
    if wtype == 'comdy':
        qids = local.day_all
        fivar = aobj.fivarvd
    
    var = 0.
    for q in qids:
        ivars = enmap.read_map( local.init_analysis_params(qid=q).fivar )
        ivars = enmap.to_healpix(ivars[0],nside=aobj.nside)
        ivars[ivars==0] = np.inf
        var += 1/ivars
    var[var==0] = np.inf
    ivar = np.nan_to_num(1./var) * amask
    print(np.min(ivar),np.max(ivar))

    hp.fitsfunc.write_map(fivar,ivar,overwrite=True)
    
