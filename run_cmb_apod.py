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
    #'snmax': 0, \
    'wind': 'base', \
    #'wind': 'com16', \
    #'ivar': 'noivar', \
    'ivar': 'base', \
    'fltr': 'none', \
    'ptsr': 'base', \
}

#//// for combined data set ////#

#qids = ['boss_d03']
#cqid = None

cqid = 'boss_s16d'
#cqid = 'boss_s15d'
#cqid = 'boss_s15n'
qids = local.get_subqids(cqid)

#run = ['alm']
run = ['alm','aps']
#run = ['aps']
#run = []

interface(run,qids,cqid=cqid,kwargs_ov=kwargs_ov,kwargs=kwargs_apod)

