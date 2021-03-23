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


def aps(qobj,rlzs,fpalm,wn,verbose=True,mean_sub=True):
    # Compute aps of reconstructed lensing map

    for q in tqdm.tqdm(qobj.qlist,ncols=100,desc='aps'):

        cl = np.zeros((len(rlzs),4,qobj.olmax+1))

        W2, W4 = wn[2], wn[4]
        
        #if mean_sub:
        #    mfg, mfc = pickle.load(open(qobj.f[q].MFalm,"rb"))
        #else:
        #    mfg, mfc = 0., 0.

        for ii, rlz in enumerate(tqdm.tqdm(rlzs,ncols=100,desc='each rlz ('+q+'):')):

            # load reconstructed kappa and curl alms with mean-field correction
            glm, clm = quad_func.load_rec_alm(qobj,q,rlz,mean_sub=mean_sub)
            #glm, clm = pickle.load(open(qobj.f[q].alm[rlz],"rb"))
            #if meansub:
            #    mfg, mfc = pickle.load(open(qobj.f[q].mfalm[rlz],"rb"))
            #else:
            #    mfg, mfc = 0., 0.

            # load kappa
            if rlz != 0:
                klm = load_input_plm(fpalm[rlz],qobj.olmax)
            else:
                klm = 0.*glm
            
            # compute cls
            cl[ii,0,:] = curvedsky.utils.alm2cl(qobj.olmax,glm)/W4
            cl[ii,1,:] = curvedsky.utils.alm2cl(qobj.olmax,clm)/W4
            cl[ii,2,:] = curvedsky.utils.alm2cl(qobj.olmax,glm,klm)/W2
            cl[ii,3,:] = curvedsky.utils.alm2cl(qobj.olmax,klm)
            np.savetxt(qobj.f[q].cl[rlz],np.concatenate((qobj.l[None,:],cl[ii,:,:])).T)

        # save sim mean
        if rlzs[1]>=1 and len(rlzs)>1:
            np.savetxt(qobj.f[q].mcls,np.concatenate((qobj.l[None,:],np.mean(cl[1:,:,:],axis=0),np.std(cl[1:,:,:],axis=0))).T)



def interface(qid,run=['norm','qrec','n0','mean','aps'],mean_sub=True,kwargs_ov={},kwargs_cmb={},kwargs_qrec={}):

    aobj = local.init_analysis_params(qid=qid,**kwargs_cmb)
    
    if qid == 'diff_dn':
        aobj_c = local.init_analysis_params(qid='comb_dn',**kwargs_cmb)
    else:
        # same as aobj
        aobj_c = local.init_analysis_params(qid=qid,**kwargs_cmb)

    # load wfactors
    wn = tools_cmb.get_wfactors([aobj.qid],aobj.ascale,wind=aobj.wind,ivar=aobj.ivar,ptsr=aobj.ptsr,fltr=aobj.fltr)[aobj.qid]

    if aobj.fltr == 'cinv':

        # filter
        ep  = 1e-30
        
        # noise
        if aobj.qid in local.wqids:
            nl = 0.
            qids = local.get_subqids(aobj.qid)
            for q in qids:
                if q in local.boss_d:
                    bobj = local.init_analysis_params(qid=q,fltr='none',wind='com16',ivar='base')
                if q in local.boss_n or q in local.s_16_d:
                    bobj = local.init_analysis_params(qid=q,fltr='none',wind='base',ivar='base')
                nl += 1. / ( np.loadtxt(bobj.fscl['n'],unpack=True)[1] + ep )
            nl = 1./(nl+ep)
        else:
            # white
            bl  = tools_cmb.beam_func(aobj.lmax,aobj.qid)
            nl  = ( local.qid_wnoise(qid) / bl )**2
        
        # corrected factors
        cnl = aobj.lcl[0,:] + nl
        #wcl = (np.loadtxt(aobj.fscl['c'])).T[1] # wiener-fileterd CMB aps
        #ocl, ifl = quad_func.cinv_empirical_fltr(aobj.lcl[0,:],wcl,cnl)
        #ocl = np.reshape( aobj.lcl[0,:]**2/(wcl+ep) ,(1,aobj.lmax+1) )
        ocl = np.reshape( cnl ,(1,aobj.lmax+1) )
        ifl = np.reshape( aobj.lcl[0,:], (1,aobj.lmax+1) )

        wn[:] = wn[0]

    else:

        # filter
        ocl = np.ones((3,aobj_c.lmax+1))
        ocl[0,:] = (np.loadtxt(aobj_c.fscl['c'])).T[1]
        ifl = ocl.copy()
    

    dirs = local.data_directory()
    qobj = quad_func.reconstruction(dirs['local'],aobj.ids,rlz=aobj.rlz,stag=aobj.stag,run=run,wn=wn,lcl=aobj.lcl,ocl=ocl,ifl=ifl,falm=aobj.falm['c'],**kwargs_ov,**kwargs_qrec)

    # Aps of reconstructed phi
    if 'aps' in run:
        aps(qobj,aobj.rlz,aobj.fiklm,wn,mean_sub=mean_sub)



