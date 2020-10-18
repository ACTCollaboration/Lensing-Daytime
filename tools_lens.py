# from external
import numpy as np
import healpy as hp
import sys
import pickle
import tqdm

# from cmblensplus/wrap/
import curvedsky

# from cmblensplus/utils/
import misctools
import quad_func

# local
import local
import tools_cmb


def load_input_plm(fpalm,lmax,verbose=False):

    if verbose: print('load input phi alms') # true phi

    # load input phi alms
    alm = np.complex128(hp.fitsfunc.read_alm(fpalm))
    # convert order of (l,m) to healpix format
    alm = curvedsky.utils.lm_healpy2healpix(alm,5100)[:lmax+1,:lmax+1]
    # convert to kappa
    L = np.linspace(0,lmax,lmax+1)
    alm = L[:,None]*(L[:,None]+1)*alm/2.

    return alm


def aps(qobj,rlzs,fpalm,wn,verbose=True,meansub=True):
    # Compute aps of reconstructed lensing map

    for q in tqdm.tqdm(qobj.qlist,ncols=100,desc='aps'):

        cl = np.zeros((len(rlzs),4,qobj.olmax+1))

        W2, W4 = wn[2], wn[4]

        for ii, rlz in enumerate(tqdm.tqdm(rlzs,ncols=100,desc='each rlz ('+q+'):')):

            # load reconstructed kappa and curl alms
            glm, clm = pickle.load(open(qobj.f[q].alm[rlz],"rb"))
            if meansub:
                mfg, mfc = pickle.load(open(qobj.f[q].mfb[rlz],"rb"))
            else:
                mfg, mfc = 0., 0.

            # load kappa
            if rlz != 0:
                klm = load_input_plm(fpalm[rlz],qobj.olmax)
            else:
                klm = 0.*glm

            # compute cls
            cl[ii,0,:] = curvedsky.utils.alm2cl(qobj.olmax,glm-mfg)/W4
            cl[ii,1,:] = curvedsky.utils.alm2cl(qobj.olmax,clm-mfc)/W4
            cl[ii,2,:] = curvedsky.utils.alm2cl(qobj.olmax,glm-mfg,klm)/W2
            cl[ii,3,:] = curvedsky.utils.alm2cl(qobj.olmax,klm)
            np.savetxt(qobj.f[q].cl[rlz],np.concatenate((qobj.l[None,:],cl[ii,:,:])).T)

        # save sim mean
        if rlzs[1]>=1 and len(rlzs)>1:
            np.savetxt(qobj.f[q].mcls,np.concatenate((qobj.l[None,:],np.mean(cl[1:,:,:],axis=0),np.std(cl[1:,:,:],axis=0))).T)



def interface(qid,run=['norm','qrec','n0','mean','aps'],meansub=True,kwargs_ov={},kwargs_cmb={},kwargs_qrec={}):

    #if ocl is None:
        # Compute filtering
        #ocl = np.ones((3,aobj.lmax+1))
        #ocl[0,:] = (np.loadtxt(aobj.fscl['c'])).T[1]
    
        # observed cl (no suppression)
        #Aobj = {q: local.init_analysis_params(qid=q,ascale=aobj.ascale) for q in qids}
        #ncl  = {q: (np.loadtxt(Aobj[q].fscl['n'])).T[1] for q in qids}
        #Ncl  = tools_cmb.comb_Nl(qids,ncl)
        #ocl  = Aobj.lcl[0,:] + Ncl

    aobj = local.init_analysis_params(qid=qid,**kwargs_cmb)
    
    if qid == 'diff_dn':
        aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs_cmb)
    else:
        aobj_c = local.init_analysis_params(qid=qid,**kwargs_cmb)

    ocl = np.ones((3,aobj_c.lmax+1))
    ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]

    # load wfactors
    wn = tools_cmb.get_wfactors([aobj.qid],aobj.ascale,wtype=aobj.wtype)[aobj.qid]

    ifl = ocl#p.lcl[0:3,:]

    dirs = local.data_directory()
    qobj = quad_func.reconstruction(dirs['local'],aobj.ids,rlz=aobj.rlz,stag=aobj.stag,run=run,wn=wn,lcl=aobj.lcl,ocl=ocl,ifl=ifl,falm=aobj.falm['c'],**kwargs_ov,**kwargs_qrec)

    # Aps of reconstructed phi
    if 'aps' in run:
        aps(qobj,aobj.rlz,aobj.fiklm,wn,meansub=meansub)



