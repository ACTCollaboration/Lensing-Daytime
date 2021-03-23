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


#qids = ['boss_dcom15','boss_com15','boss_dcom16','boss_com']

aobj = {q: local.init_analysis_params(qid=q,ascale=ap) for q in local.qid_all}

mask = {}

for q in qids:
    mask_iv = tools_cmb.load_mask(q,with_ivar=False)
    mask[q] = enmap.to_healpix(mask_iv,nside=2048)

if wqid == 'boss_dcom15':
    mask['boss_d03']
    hp.fitsfunc.write_map(aobj_dcom15.amask,mask,overwrite=True)

