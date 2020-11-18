#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np, curvedsky, pickle
from pixell import enmap
from tools_cmb import *



kwargs_cinv = {\
    'chn' : 1, \
    'eps' : [1e-5], \
    'lmaxs' : [4096], \
    'nsides' : [2048], \
    'itns' : [100], \
    'ro' : 1, \
    'filter' : 'W'
}

kwargs_cmb = {\
    'snmin': 0, \
    'snmax': 10
}

wiener_cinv_core('boss_d03',kwargs_cmb=kwargs_cmb,kwargs_cinv=kwargs_cinv)


