#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
from tools_cmb import *
import tools_lens
import local


kwargs_ov = {\
    'overwrite':True, \
    'verbose':True
}

kwargs_cmb = {\
    'snmin': 0, \
    'snmax': 100, \
    'fltr': 'cinv' \
}

cqid = 'boss_s15dn'

# for comparison
interface(['aps'],[cqid],kwargs_ov=kwargs_ov,kwargs=kwargs_cmb)

