
import numpy as np, os
from quicksub import *

def create_runfile(f,qid,snmin=0,snmax=100,ascale=1.,overwrite=False,verbose=False):
    
    add('import numpy as np, tools_cmb',f,ini=True)
    add("kwargs = {'snmin':"+str(snmin)+",'snmax':"+str(snmax)+",'ascale':"+str(ascale)+"}",f)
    if isinstance(qid,str):
        add("qids = ['"+qid+"']",f)
    if isinstance(qid,list):
        add("qids = ['"+qid[0]+"','"+qid[1]+"']",f)
    #add("tools_cmb.generate_map(qids,overwrite="+str(overwrite)+",verbose=True,**kwargs)",f)
    add("tools_cmb.map2alm(qids,overwrite="+str(overwrite)+",verbose="+str(verbose)+",**kwargs)",f)
    add("tools_cmb.alm2aps(qids,overwrite="+str(overwrite)+",verbose="+str(verbose)+",mtype=['T'],**kwargs)",f)
    add("tools_cmb.alm_supfac(qids,overwrite="+str(overwrite)+",verbose="+str(verbose)+",**kwargs)",f)


def jobfile(tag,qid,**kwargs):

    #set run file
    f_run = 'tmp_job_run_'+tag+'.py'
    create_runfile(f_run,qid,**kwargs)

    # set job file
    f_sub = 'tmp_job_sub_'+tag+'.sh'
    set_sbatch_params(f_sub,tag,mem='64G',t='0-12:00',email=True)
    add('source ~/.bashrc.ext',f_sub)
    add('py4so',f_sub)
    add('python '+f_run,f_sub)
    
    # submit
    #os.system('sbatch '+f_sub)
    os.system('sh '+f_sub)
    os.system('rm -rf '+f_run+' '+f_sub)


#qids = ['boss_d01','boss_d02',['boss_d03','boss_d04'],'boss_01','boss_02',['boss_03','boss_04'],'s16_d01',['s16_d02','s16_d03']]
#qids = ['boss_d01','boss_d02','boss_01','boss_02','s16_d01']
#qids = [['s16_d02','s16_d03']]
qids = [['boss_d03','boss_d04']]
for qid in qids:
    if isinstance(qid,str): tag = qid
    if isinstance(qid,list): tag = '_'.join(qid)
    jobfile(tag,qid,overwrite=True,verbose=False,snmin=0,snmax=100,ascale=5.)

