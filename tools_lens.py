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


def load_input_plm(fpalm,lmax,verbose=False):

    if verbose: print('load input phi alms') # true phi

    # load input phi alms
    alm = np.complex128(hp.fitsfunc.read_alm(fpalm))
    # convert order of (l,m) to healpix format
    alm = curvedsky.utils.lm_healpy2healpix(len(alm),alm,5100)[:lmax+1,:lmax+1]

    return alm


def aps(qobj,rlzs,fpalm,wn,verbose=True):
    # Compute aps of reconstructed lensing map

    for q in tqdm.tqdm(qobj.qlist,ncols=100,desc='aps'):

        cl = np.zeros((len(rlzs),4,qobj.olmax+1))

        W2, W4 = wn[2], wn[4]

        for ii, rlz in enumerate(tqdm.tqdm(rlzs,ncols=100,desc='each rlz ('+q+'):')):

            # load reconstructed kappa and curl alms
            glm, clm = pickle.load(open(qobj.f[q].alm[rlz],"rb"))
            mfg, mfc = pickle.load(open(qobj.f[q].mfb[rlz],"rb"))

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


def interface(aobj,wn,run=['norm','qrec','n0','mean','aps'],ocl=None,kwargs_ov={},kwargs_qrec={}):

    if ocl is None
        # Compute filtering
        ocl = np.ones((3,aobj.lmax+1))
        ocl[0,:] = (np.loadtxt(aobj.fscl['c'])).T[1]
    
    ifl = ocl#p.lcl[0:3,:]

    dirs = local.data_directory()
    qobj = quad_func.reconstruction(dirs['local'],aobj.ids,rlz=aobj.rlz,stag=aobj.stag,run=run,wn=wn,lcl=aobj.lcl,ocl=ocl,ifl=ifl,falm=aobj.falm['c'],**kwargs_ov,**kwargs_qrec)

    # Aps of reconstructed phi
    if 'aps' in run:
        aps(qobj,aobj.rlz,aobj.fiklm,wn)



