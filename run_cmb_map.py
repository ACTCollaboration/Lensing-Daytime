#!/usr/bin/env python
# coding: utf-8

import numpy as np
import local
import tools_cmb

# global fixed parameters
ow = True
vb = False

qids = local.qid_all
#qids = ['boss_d01','boss_d02','boss_d03']
kwargs = {'snmin':0,'snmax':0}

tools_cmb.generate_map(qids,overwrite=ow,verbose=vb,**kwargs)

