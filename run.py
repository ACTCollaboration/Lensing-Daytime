#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import tools_cmbalm


# In[2]:

ow = True

qids = ['boss_d01','boss_d02','boss_d03','boss_d04','boss_01','boss_02','boss_03','boss_04','s16_d01','s16_d02','s16_d03']
#qids = ['s16_d02','s16_d03']
#qids = ['boss_04']

kwargs = {'snmin':0,'snmax':10}

print('Generating map')
tools_cmbalm.generate_map(qids,overwrite=ow,verbose=True,**kwargs)

print('Convert map to alm')
tools_cmbalm.map2alm(qids,overwrite=ow,**kwargs)

print('Compute aps')
tools_cmbalm.alm2aps(qids,overwrite=ow,verbose=True,mtype=['T'],**kwargs)

print('Compute dif and cross spectra')
for q in qids:
    tools_cmbalm.alm2aps_null(q,overwrite=ow,verbose=True,**kwargs)

