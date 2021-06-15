#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
from tools_cmb import *
import tools_lens
import local


kwargs_ov = {\
    'overwrite':False, \
    'verbose':True
}

kwargs_cinv = {\
    'chn' : 1, \
    'eps' : [1e-4], \
    'lmaxs' : [4096], \
    'nsides' : [2048], \
    'itns' : [500], \
    'ro' : 1, \
    'filter' : 'W', \
    'stat': 'status1.txt' \
}

kwargs_cmb = {\
    'snmin': 0, \
    'snmax': 100, \
    'fltr': 'cinv', \
    'ptsr': 'base', \
}

wqid = 'boss_s15n'
qids = local.get_subqids(wqid)

# cinv
wiener_cinv_core(qids,wqid,kwargs_ov=kwargs_ov,kwargs_cmb=kwargs_cmb,kwargs_cinv=kwargs_cinv)
interface(['aps'],[wqid],kwargs_ov=kwargs_ov,kwargs=kwargs_cmb)

