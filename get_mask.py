#!/usr/bin/env python
# coding: utf-8

import numpy as np, healpy as hp, curvedsky as cs, sys, local, tools_cmb

sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/actsims/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/soapack/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/orphics/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/tilec/")
from pixell import enmap
from soapack import interfaces

#run = ['apod']
#run = ['iso']
#run = ['com']
run = ['PTSR']
ptsr = 'base'
#run = ['cinv']

    
if 'apod' in run:

    #qids = ['boss_d01']
    qids = ['boss_01','boss_d02','boss_02','boss_d03','boss_03','boss_d04','boss_04','s16_d01','s16_d02','s16_d03']
    
    # for different apodization scales
    for ap in [1.,2.,3.,4.,5.]:
    #or ap in [1.,5.]:

        aobj = {q: local.init_analysis_params(qid=q,ascale=ap) for q in qids}

        for q in qids:
        
            mask_iv = tools_cmb.load_survey_mask(q)
            mask = enmap.to_healpix(mask_iv[0],nside=2048)
    
            if ap != 1.:
                mask[mask!=0] = 1. # binary mask
                amask = cs.utils.apodize(2048,mask,aobj[q].ascale)
            else:
                amask = mask.copy()
        
            print(q,np.average(amask),np.max(amask))        
            hp.fitsfunc.write_map(aobj[q].amask,amask,overwrite=True)


if 'cinv' in run:

    # for different apodization scales
    for wqid in local.wqids:

        aobj = local.init_analysis_params(qid=wqid,ascale=0.)
        window = 1.
        
        for q in local.get_subqids(wqid):
        
            mask = tools_cmb.load_survey_mask(q)
            Mask = enmap.to_healpix(mask,nside=2048)
            Mask[Mask!=0.] = 1.
            print(q,np.average(Mask),np.max(Mask))
            window *= (1.-Mask)

        # binary
        window = 1.-window
        print('save to '+aobj.amask)
        print(wqid,np.average(window),np.max(window))
        hp.fitsfunc.write_map(aobj.amask,window,overwrite=True)

        # apodization
        aobj = local.init_analysis_params(qid=wqid,ascale=1.)
        amask  = cs.utils.apodize(2048,window,aobj.ascale)
        print('save to '+aobj.amask)
        hp.fitsfunc.write_map(aobj.amask,amask,overwrite=True)

            
if 'com' in run: # overlapped region

    for wtype in ['com15','com16']:
        
        if wtype == 'com15':
            qids = local.boss_dn
        if wtype == 'com16':
            qids = local.qid_all

        # compute common binary region
        common = 1.
        for q in qids:
            mask = tools_cmb.load_survey_mask(q)
            mask[mask!=0] = 1.
            common *= enmap.to_healpix(mask,nside=2048)

        #for ap in [3.,5.]:
        for ap in [1.]:
            amask = cs.utils.apodize(2048,common,ap)
            aobj = local.init_analysis_params(ascale=ap,wtype=wtype)
            hp.fitsfunc.write_map(aobj.amask,amask,overwrite=True)


            
if 'iso' in run: 

    qids_s15 = local.boss_dn
    qids_s16 = local.s_16_d

    # compute s15 only binary region
    common = {}
    for qid, qids in [('s15',qids_s15),('s16',qids_s16)]:
        common[qid] = 1.
        for q in qids:
            mask = tools_cmb.load_survey_mask(q)
            mask[mask!=0] = 1.
            common[qid] *= enmap.to_healpix(mask,nside=2048)
            
    isomask = common['s15'] * (1-common['s16'])
    #for ap in [3.,5.]:
    for ap in [1.]:
        amask = cs.utils.apodize(2048,isomask,ap)
        aobj = local.init_analysis_params(ascale=ap,wtype='iso15')
        if ap==.5:
            amask[amask<.5] = 0.
            amask[amask!=0] = 1.
            amask = cs.utils.apodize(2048,amask,3.)
        hp.fitsfunc.write_map(aobj.amask,amask,overwrite=True)


#if 'ptsr' in run:
    
#    print('creating old ptsr mask')
    
#    aobj = local.init_analysis_params()
#    ptsr = tools_cmb.create_ptsr_mask_old(aobj.nside)
#    hp.fitsfunc.write_map(aobj.fptsr_old,ptsr,overwrite=True)


if 'PTSR' in run:
    
    print('creating ptsr mask')
    
    aobj = local.init_analysis_params(ptsr=ptsr)
    mask = tools_cmb.create_ptsr_mask(aobj.nside,ptsr=aobj.ptsr)
    hp.fitsfunc.write_map(aobj.fptsr,mask,overwrite=True)


