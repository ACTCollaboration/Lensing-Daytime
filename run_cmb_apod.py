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

kwargs_apod = {\
    'snmin': 0, \
    'snmax': 100, \
    'ivar': 'noivar', \
    'wind': 'base', \
}

#//// for combined data set ////#

#cqid = 'boss_s15d'
#cqid = 'boss_s15n'
#cqid = 'boss_s16d'
#cqid = 'boss_s15dn'
cqid = 'boss_alldn'

#run = ['alm','aps','comb']
#run = ['aps','comb']
run = ['comb']

qids = local.get_subqids(cqid)
interface(run,qids,cqid=cqid,kwargs_ov=kwargs_ov,kwargs=kwargs_apod)

