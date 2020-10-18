
# from external
import numpy as np
import healpy as hp
import sys
import pickle
import tqdm

# from act library
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/actsims/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/soapack/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/orphics/")
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/tilec/")
from actsims import simgen, noise
from pixell import enmap
from soapack import interfaces

# from cmblensplus/wrap/
import curvedsky

# from cmblensplus/utils/
import misctools

# from this directory
import local


def interface(run,qids,kwargs_ov={},kwargs={}):
    
    if 'alm' in run:
        map2alm(qids,**kwargs_ov,**kwargs)

    if 'aps' in run:
        alm2aps(qids,mtype=['T'],**kwargs_ov,**kwargs)

    if 'sup' in run:
        alm_supfac(qids,**kwargs_ov,**kwargs)
        
    if 'each' in run:
        for qid_d, qid_n in zip(local.boss_d,local.boss_n):
            diff_day_night(qid_d,qid_n,'diff_'+qid_n,mtype=['T'],**kwargs_ov,**kwargs)

    if 'comb' in run:
    
        qid_dn = ['comb_d','comb_n']
        if 'iso15' in kwargs['wtype']:
            if 'V3' in kwargs['wtype']:
                qids_d = ['boss_d02','boss_d03','boss_d04']
                qids_n = ['boss_02','boss_03','boss_04']
                #qids_d = local.boss_d
                #qids_n = local.boss_n
            else:
                qids_d = ['boss_d03']
                qids_n = ['boss_03','boss_04']
        elif 'com15' in kwargs['wtype']:
            qids_d = local.boss_d
            qids_n = local.boss_n
        else:
            qids_d = local.day_all
            qids_n = local.boss_n

        # comb day or night
        print('combine alm for day')
        alm_comb(qids_d,qid_dn[0],mtype=['T'],**kwargs_ov,**kwargs)
        print('combine alm for night')
        alm_comb(qids_n,qid_dn[1],mtype=['T'],**kwargs_ov,**kwargs)
    
        alm2aps(qid_dn,mtype=['T'],**kwargs_ov,**kwargs)
        alm_supfac(qid_dn,**kwargs_ov,**kwargs)

        # day - night
        print('compute day - night')
        diff_day_night(qid_dn[0],qid_dn[1],'diff_dn',mtype=['T'],**kwargs_ov,**kwargs)

        # day + night
        alm_comb(qids_d+qids_n,'comb_dn',mtype=['T'],**kwargs_ov,**kwargs)
        alm2aps(['comb_dn'],mtype=['T'],**kwargs_ov,**kwargs)


def coadd_real_data(qid,mask,dg=2):
    # Compute coadded map by simply summing up each split with the inverse variance weight 

    assert mask.ndim == 2

    # First, read survey parameters from qid
    model, season, array, patch, freq = local.qid_info(qid[0])
    
    # Day time data
    if model == 'dr5':
        dm = interfaces.models['dr5'](region=mask)
        datas = enmap.enmap( [ dm.get_splits( q, calibrated=True ) for q in qid ] )
        ivars = enmap.enmap( [ dm.get_ivars( q, calibrated=True ) for q in qid ] )

    # Night time data
    if model == 'act_mr3':
        dm = interfaces.models['act_mr3'](region=mask,calibrated=True)
        datas = dm.get_splits(season=season,patch=patch,arrays=dm.array_freqs[array],srcfree=True)
        ivars = dm.get_splits_ivar(season=season,patch=patch,arrays=dm.array_freqs[array])

    # normalize
    #ivars = ivars_normalize(ivars)
    
    # coadded map
    civars = np.average(ivars,axis=1) 
    civars[civars==0] = np.inf
    map_c = mask[None,None,:,:] * np.average( datas*ivars, axis=1 )/ civars / local.Tcmb 
    
    # downgrade
    map_c = enmap.downgrade(map_c, dg)
    
    # check
    assert map_c.ndim == 4

    return map_c


def get_ivar(qid,mask,dg=2):

    # First, read survey parameters from qid
    model, season, array, patch, freq = local.qid_info(qid[0])

    # Day time data
    if model == 'dr5':
        dm = interfaces.models['dr5'](region=mask)
        ivars = enmap.enmap( [ dm.get_ivars( q, calibrated=True ) for q in qid ] )

    # Night time data
    if model == 'act_mr3':
        dm = interfaces.models['act_mr3'](region=mask,calibrated=True)
        ivars = dm.get_splits_ivar(season=season,patch=patch,arrays=dm.array_freqs[array])

    # downgrade
    ivars = enmap.downgrade(ivars, dg)

    return ivars


def ivars_normalize(ivars):
    
    n1, n2, n3 = np.shape(ivars[:,:,:,0,0])
    nivars = 0.*ivars
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                nivars[i,j,k,:,:] = ivars[i,j,k,:,:]/np.max(ivars[i,j,k,:,:])
    return nivars
    

def generate_map(qids,overwrite=False,verbose=True,dg=2,**kwargs):
    # Here we compute the real coadd map and simulated maps from the noise covariance and pre-computed fullsky signal
    # This function uses actsims "simgen" code

    for qid in qids:
        
        if qid in ['boss_d04','s16_d03','boss_04']:
            # Here, the simulation is performed for each array, not frequency, to take into account their correlation
            # For pa3, there are two frequency bands
            if verbose:
                print('skip '+qid+': the map is generated by running f090 case.')
            continue
        
        if '_d0' in qid:
            version = 'v6.3.0_calibrated_mask_version_masks_20200723'
        else:
            version = 'v6.3.0_calibrated_mask_version_padded_v1'

        # define qid_array to take into account multi-frequency case 
        if qid == 'boss_d03':
            qid_array = ['boss_d03','boss_d04']
        elif qid == 's16_d02':
            qid_array = ['s16_d02','s16_d03']
        elif qid == 'boss_03':
            qid_array = ['boss_03','boss_04']            
        else:
            # single frequency case
            qid_array = [qid]

        # define filename
        aobj = { q: local.init_analysis_params(qid=q,**kwargs) for q in qid_array }

        # load survey mask
        mask = load_survey_mask(qid,dg=1)[0]
        
        # Define an object for sim generation
        model, season, array, patch, __ = local.qid_info(qid)
        if model == 'dr5':
            simobj = simgen.SimGen(version=version,qid=qid_array,model=model)
        if model == 'act_mr3':
            simobj = simgen.SimGen(version=version,model=model)

        
        # save 1/var map
        if not misctools.check_path(aobj[qid_array[0]].fivar,overwrite=overwrite,verbose=verbose): 
            ivars = get_ivar(qid_array,mask)
            civars = np.average(ivars,axis=1) # coadd ivar
            for qi, q in enumerate(qid_array):
                enmap.write_map( aobj[q].fivar, civars[qi,:,:,:] )


        # loop over realizations
        for i in tqdm.tqdm(aobj[qid].rlz):
            
            if misctools.check_path(aobj[qid].fmap['s'][i],overwrite=overwrite,verbose=verbose): continue

            # maps will have signal ('s') and noise ('n') and each has [Nfreq, NTPol, Coord1, Coord2]
            maps = {}

            # real data
            if i == 0: 
                maps['s'] = coadd_real_data(qid_array,mask) 
            # simulation
            else:
                maps['s'], maps['n'], ivars = simobj.get_sim(season, patch, array, sim_num=i, set_idx=0)

                # normalize
                #ivars = ivars_normalize(ivars)

                # coadd with civars weight implicitly multiplied
                civars = np.average(ivars,axis=1)
                civars[civars==0] = np.inf
                maps['s'] = mask[None,None,:,:] * np.average( maps['s']*ivars, axis=1 ) / civars / local.Tcmb 
                maps['n'] = mask[None,None,:,:] * np.average( maps['n']*ivars, axis=1 ) / civars / local.Tcmb

                # downgrade
                maps['s'] = enmap.downgrade(maps['s'], dg)
                maps['n'] = enmap.downgrade(maps['n'], dg)

            # save signal and noise to files
            # data is saved to a signal file
            for s in ['s','n']:
                if s=='n' and i==0: continue
                for qi, q in enumerate(qid_array):
                    enmap.write_map(aobj[q].fmap[s][i],maps[s][qi,:,:,:])


                
def load_survey_mask_core(qid):    
    # load mask by specifying qid
    
    model, season, array, patch, freq = local.qid_info(qid)
    
    if model == 'dr5':
        mask = interfaces.get_binary_apodized_mask(qid)
    
    elif model == 'act_mr3':
        mask = interfaces.get_act_mr3_crosslinked_mask(patch,version='padded_v1',kind='binary_apod',season=season,array=array+"_"+freq,pad=None)

    return mask


def load_survey_mask(qid,dg=2):

    # load survey mask
    mask = load_survey_mask_core(qid)
    mask_dg = enmap.downgrade(mask,dg)

    return mask_dg


def load_ivar_curvedsky(aobj):
    
    wtype = (aobj.wtype).replace('pt','')
    
    if 'vc' in wtype:
        if '16' in wtype:
            return hp.fitsfunc.read_map(aobj.fivar16,verbose=False)
        elif '15' in wtype:
            return hp.fitsfunc.read_map(aobj.fivar15,verbose=False)
        else:
            sys.exit('ivar is not well defined')
    elif 'vd' in wtype:
        ivar = enmap.read_map( aobj.fivarvd )
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    elif 'v3' in wtype:
        ivar = enmap.read_map( local.init_analysis_params(qid='boss_d03').fivar )
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    elif 'V3' in wtype:
        ivar = enmap.read_map( local.init_analysis_params(qid='boss_d03').fivar ) ** 2
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    elif 'v0' in wtype:
        return 1.
    elif wtype == 'base':
        ivar = enmap.read_map( aobj.fivar )
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    else:
        sys.exit('ivar is not well defined')
        

def create_square_mask(nside,lon,lat,theta,phi):
    mask = np.zeros(12*nside**2)
    mask[lon[0]>=theta] = 1.
    mask[lon[1]<=theta] = 1.
    mask[lat[0]>=phi] = 1.
    mask[lat[1]<=phi] = 1.
    return mask


def create_ptsr_mask_old(nside,\
                lonras = [[186.8,187.8],[187.2,188.2],[208.8,209.8],[237,237.7],[164.3,164.9],[225.7,226.3]],\
                latras = [[1.5,2.5],[12.1,12.8],[19,19.7],[2.3,3],[1.3,1.9],[10.2,10.8]],\
                #lonras = [[186.8,187.8],[187.2,188.2],[164.3,164.9],[225.7,226.3]],\
                #latras = [[1.5,2.5],[12.1,12.8],[1.3,1.9],[10.2,10.8]],\
               ascale=0.5):
    
    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(12*nside**2), lonlat=True)
    maskpt = 1.
    
    for lonra, latra in zip(lonras,latras):
        maskpt *= create_square_mask(nside,lonra,latra,pixel_theta,pixel_phi)
    if ascale!= 0.:
        maskpt = curvedsky.utils.apodize(nside,maskpt,ascale)

    return maskpt


def create_ptsr_mask(nside,ascale=0.5,threshold=.0,extend=3.):
    
    # from catalogue
    #ras, decs, size = np.loadtxt('data_local/input/cat_crossmatched.txt',unpack=True,usecols=(0,1,5))
    #RAs  = ras[size>threshold]
    #DECs = decs[size>threshold]
    #arcm = size[size>threshold] * extend
    
    RAs, DECs = interfaces.get_act_mr3f_union_sources(version='20200503_sncut_40')
    arcm = np.ones(len(RAs)) * 10. * extend

    # additional mask
    RAs_add  = np.array([187.3])
    DECs_add = np.array([2.])
    arcm_add = np.array([40.])
    #RAs_add  = np.array([187.3,187.7,209.3,164.6,226.])
    #DECs_add = np.array([2.,12.5,19.3,1.6,10.5])
    #arcm_add = np.array([30.,30.,10.,10.,10.])
    # add
    RAs  = np.concatenate((RAs,RAs_add))
    DECs = np.concatenate((DECs,DECs_add))
    arcm = np.concatenate((arcm,arcm_add))
        
    # compute 3D positions
    v = hp.pixelfunc.ang2vec(RAs, DECs, lonlat=True)
    
    # create mask
    maskpt = np.ones(12*nside**2)
    
    for i in range(len(arcm)):
        pix = hp.query_disc(nside, v[i], arcm[i]*np.pi/10800.)
        maskpt[pix] = 0.

    if ascale!= 0.:
        maskpt = curvedsky.utils.apodize(nside,maskpt,ascale)

    return maskpt


def load_window_curvedsky( aobj, survey_mask=True, with_ivar=True, ivar_norm=False ):

    # load survey mask
    if survey_mask and ( aobj.qid in local.qid_all ):
        mask = load_survey_mask(aobj.qid)
        mask = enmap.to_healpix(mask,nside=aobj.nside)  # to healpix map
    else:
        mask = 1.

    # load ivar
    if with_ivar:
        ivar = load_ivar_curvedsky(aobj)
        if ivar_norm:
            ivar = ivar/np.max(ivar)
    else:
        ivar = 1.

    # load additional apodization mask
    if ( aobj.ascale!=1. and aobj.qid in local.qid_all ) or 'com' in aobj.wtype or 'iso' in aobj.wtype:
        amask = hp.fitsfunc.read_map( aobj.amask, verbose=False )
    #elif aobj.ascale!=1. and aobj.qid not in local.qid_all:
    #    amask = hp.fitsfunc.read_map( local.init_analysis_params(qid='boss_03').amask ,verbose=False)
    else:
        amask = 1.

    # load ptsr mask for day
    if 'pt' in aobj.wtype:
        ptsr = hp.fitsfunc.read_map( aobj.fptsr_old ,verbose=False )
        #ptsr = create_ptsr_mask_old(aobj.nside)
    elif 'PT' in aobj.wtype:
        ptsr = hp.fitsfunc.read_map( aobj.fptsr ,verbose=False )
        #ptsr = create_ptsr_mask(aobj.nside)
    else:
        ptsr = 1.

    return mask*ivar*amask*ptsr

 
def get_wfactor(mask,wnmax=5):
    wn = np.zeros(wnmax)
    for n in range(1,wnmax):
        wn[n] = np.average(mask**n)
    wn[0] = np.average(mask/(mask+1e-30))
    print('wfactors:',wn)
    return wn
    

def get_wfactors(qids,ascale,wtype='base',with_ivar=True,ivar_norm=False,wnmax=5):
    w = {}
    for q in qids:
        aobj = local.init_analysis_params(qid=q,ascale=ascale,wtype=wtype)
        mask = load_window_curvedsky(aobj,with_ivar=with_ivar,ivar_norm=ivar_norm)
        w[q] = get_wfactor(mask,wnmax=wnmax)
    return w

   
def define_lmask(shape,wcs, lxcut = None, lycut = None, lmin = None, lmax = None):
    # copied from orphics
    output = np.ones(shape[-2:], dtype = int)
    if (lmin is not None) or (lmax is not None): modlmap = enmap.modlmap(shape, wcs)
    if (lxcut is not None) or (lycut is not None): ly, lx = enmap.laxes(shape, wcs, oversample=1)
    if lmin is not None:
        output[np.where(modlmap <= lmin)] = 0
    if lmax is not None:
        output[np.where(modlmap >= lmax)] = 0
    if lxcut is not None:
        output[:,np.where(np.abs(lx) < lxcut)] = 0
    if lycut is not None:
        output[np.where(np.abs(ly) < lycut),:] = 0
    return output
    

def remove_lxly(fmap,lmin=100,lmax=4096):
    
    alm = enmap.fft(fmap)
    shape, wcs = fmap.shape, fmap.wcs
    kmask = define_lmask(shape,wcs,lmin=lmin,lmax=lmax,lxcut=90,lycut=50)
    #print(np.shape(kmask),np.shape(alm))
    alm[kmask<0.5] = 0
    fmap_fl = enmap.ifft(alm).real
    return fmap_fl


def map2alm_core(fmap,aobj,amask):

    fmap_fl = remove_lxly( fmap, lmin=aobj.clmin, lmax=aobj.lmax )
    hpmap = enmap.to_healpix( fmap_fl, nside=aobj.nside )
    alm = curvedsky.utils.hp_map2alm( aobj.nside, aobj.lmax, aobj.lmax, hpmap*amask )

    return alm


def beam_func(lmax,qid):
    
    dm = interfaces.models['dr5']()
    L  = np.arange(0, lmax+100, dtype=np.float)
    Bl = dm.get_beam_func(L, qid)
    return Bl[:lmax+1]


def map2alm(qids,overwrite=False,verbose=True,ascale=1.,**kwargs):

    
    # loop over qid
    for q in qids:
        
        aobj = local.init_analysis_params(qid=q,ascale=ascale,**kwargs)
        
        # beam
        Bl = beam_func(aobj.lmax,q)

        # total mask
        mask = load_window_curvedsky( aobj, survey_mask=False )
        #civar = load_ivar_curvedsky(aobj)

        # change apodization
        #if ascale != 1. or 'com' in aobj.wtype or 'iso' in aobj.wtype:
        #    if verbose: print('correct apodization mask')
        #    amask = hp.fitsfunc.read_map(aobj.amask,verbose=False)
            
        # load ptsr mask for day
        #if 'pt' in aobj.wtype:
        #    ptsr = create_ptsr_mask(aobj.nside)
        #else:
        #    ptsr = 1.

        # loop for map -> alm
        for i in tqdm.tqdm(aobj.rlz,desc='map -> alm'):
        
            if misctools.check_path(aobj.falm['c']['T'][i],overwrite=overwrite,verbose=verbose): continue
            
            if i == 0: 

                maps_c = enmap.read_map(aobj.fmap['s'][i])[0]
                alm_c = map2alm_core( maps_c, aobj, mask ) / Bl[:,None]
                pickle.dump((alm_c),open(aobj.falm['c']['T'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

            else:
                
                maps_s = enmap.read_map(aobj.fmap['s'][i])[0]
                maps_n = enmap.read_map(aobj.fmap['n'][i])[0]

                alm_s = map2alm_core( maps_s, aobj, mask ) / Bl[:,None]
                alm_n = map2alm_core( maps_n, aobj, mask ) / Bl[:,None]

                pickle.dump((alm_n),open(aobj.falm['n']['T'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump((alm_s+alm_n),open(aobj.falm['c']['T'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)



def alm2aps_core(lmax,falm,w2=1.,mtype=['T']):

    cl = np.zeros((6,lmax+1))
    
    # compute cls
    if 'T' in mtype:  
        Talm = pickle.load(open(falm['T'],"rb"))
        cl[0,:] = curvedsky.utils.alm2cl(lmax,Talm)

    if 'E' in mtype:  
        Ealm = pickle.load(open(falm['E'],"rb"))
        cl[1,:] = curvedsky.utils.alm2cl(lmax,Ealm)

    if 'B' in mtype:  
        Balm = pickle.load(open(falm['B'],"rb"))
        cl[2,:] = curvedsky.utils.alm2cl(lmax,Balm)

    if 'T' in mtype and 'E' in mtype:
        cl[3,:] = curvedsky.utils.alm2cl(lmax,Talm,Ealm)

    if 'T' in mtype and 'B' in mtype:
        cl[4,:] = curvedsky.utils.alm2cl(lmax,Talm,Balm)

    if 'E' in mtype and 'B' in mtype:
        cl[5,:] = curvedsky.utils.alm2cl(lmax,Ealm,Balm)

    return cl/w2
        

def alm2aps(qids,overwrite=False,verbose=True,mtype=['T'],W2=None,**kwargs):

    for qid in qids: 
        
        aobj = local.init_analysis_params(qid=qid,**kwargs)
        
        if W2 is None:
            mask = load_window_curvedsky(aobj)
            w2 = get_wfactor(mask)[2]
        else:
            if isinstance(W2,dict):
                w2 = W2[qid]
            else:
                w2 = W2

        cl = {s: np.zeros((len(aobj.rlz),6,aobj.lmax+1)) for s in ['c','n','s']}

        for ii, rlz in enumerate(tqdm.tqdm(aobj.rlz,desc='alm -> aps')):

            if misctools.check_path(aobj.fcls['c'][rlz],overwrite=overwrite,verbose=verbose): continue

            if rlz == 0:
                fnames = { m: aobj.falm['c'][m][rlz] for m in mtype }
                cl['c'][ii,:,:] = alm2aps_core(aobj.lmax,fnames,w2=w2,mtype=mtype)
            else:
                for s in ['c','n']:
                    fnames = { m: aobj.falm[s][m][rlz] for m in mtype }
                    cl[s][ii,:,:] = alm2aps_core(aobj.lmax,fnames,w2=w2,mtype=mtype)

                # signal part
                Tclm = pickle.load(open(aobj.falm['c']['T'][rlz],"rb"))
                Tnlm = pickle.load(open(aobj.falm['n']['T'][rlz],"rb"))
                cl['s'][ii,0,:] = curvedsky.utils.alm2cl(aobj.lmax,Tclm-Tnlm)
 
            # save cl for each rlz
            np.savetxt(aobj.fcls['c'][rlz],np.concatenate((aobj.l[None,:],cl['c'][ii,:,:])).T)
            if rlz >= 1:
                np.savetxt(aobj.fcls['s'][rlz],np.concatenate((aobj.l[None,:],cl['s'][ii,:,:])).T)
                np.savetxt(aobj.fcls['n'][rlz],np.concatenate((aobj.l[None,:],cl['n'][ii,:,:])).T)

        # save mean cl to files
        if aobj.rlz[-1] >= 2:
            if verbose:  print('save averaged spectrum over rlz')
            imin = max(0,1-aobj.rlz[0]) 
            for s in ['c','n','s']:
                if misctools.check_path(aobj.fscl[s],overwrite=overwrite,verbose=verbose): continue
                np.savetxt(aobj.fscl[s],np.concatenate((aobj.l[None,:],np.mean(cl[s][imin:,:,:],axis=0),np.std(cl[s][imin:,:,:],axis=0))).T)


def alm_supfac(qids,W1=None,overwrite=False,verbose=True,**kwargs):

    del kwargs['snmin'] # set snmin to 1 below
    
    for qid in qids:

        aobj = local.init_analysis_params(qid=qid,snmin=1,**kwargs)

        if misctools.check_path(aobj.fsup,overwrite=overwrite,verbose=verbose): continue

        if W1 is None:
            mask = load_window_curvedsky(aobj)
            w1 = get_wfactor(mask)[1]
        else:
            if isinstance(W1,dict):
                w1 = W1[qid]
            else:
                w1 = W1
            
        lmax = aobj.lmax
        rl = np.zeros((len(aobj.rlz),lmax+1))

        for ii, rlz in enumerate(tqdm.tqdm(aobj.rlz,desc='compute supfac')):
 
            # projected signal alm
            Tlm_c = pickle.load(open(aobj.falm['c']['T'][rlz],"rb"))
            Tlm_n = pickle.load(open(aobj.falm['n']['T'][rlz],"rb"))
            Tlm_obs = Tlm_c - Tlm_n
            
            # input fullsky signal alm
            f = '/project/projectdirs/act/data/actsims_data/signal_v0.4/fullskyLensedUnabberatedCMB_alm_set00_'+str(rlz).zfill(5)+'.fits'
            Tlm_inp = np.complex128( hp.fitsfunc.read_alm( f, hdu = (1) ) ) / local.Tcmb
            ilmax = hp.sphtfunc.Alm.getlmax(len(Tlm_inp))
            Tlm_inp = curvedsky.utils.lm_healpy2healpix(Tlm_inp, ilmax)[:lmax+1,:lmax+1]

            # obs x input
            xl = curvedsky.utils.alm2cl(lmax,Tlm_obs,Tlm_inp)
            
            # input auto
            cl = curvedsky.utils.alm2cl(lmax,Tlm_inp)
            
            # take ratio to get suppression factor
            rl[ii,2:] = xl[2:]/cl[2:]

        # save to file
        mrl = np.mean(rl,axis=0)
        np.savetxt( aobj.fsup, np.array( ( aobj.l, mrl, mrl/w1 ) ).T )
        


def comb_Nl(qids,ncl,rcl=None):

    Nl = 0.*ncl[qids[0]]
    
    for q in qids:
        if rcl is None:
            Nl[2:] += 1./ncl[q][2:]
        else:
            Nl[2:] += rcl[q][2:]**2/ncl[q][2:]
    Nl[2:] = 1./Nl[2:]
    return Nl



def alm_comb(qids,qidc,overwrite=False,verbose=True,mtype=['T'],ep=1e-30,**kwargs):

    # qids = qid to be combined
    # qidc = output qid
    
    aobj = {q: local.init_analysis_params(qid=q,**kwargs) for q in qids+[qidc]}

    # pre-computed spectra
    #mcl  = {q: (np.loadtxt(aobj[q].fscl['c'])).T[1] for q in qids}
    ncl  = {q: (np.loadtxt(aobj[q].fscl['n'])).T[1] for q in qids}
    # sup fac
    #fcl  = {q: (np.loadtxt(aobj[q].fsup)).T[1] for q in qids}
    #rcl  = {q: (np.loadtxt(aobj[q].fsup)).T[2] for q in qids}
    #Ncl  = comb_Nl(qids,ncl,rcl)
    Ncl  = comb_Nl(qids,ncl)

    #w1 = {}
    #for qid in qids:
    #    mask_iv = load_mask(qid)
    #    mask_hp = enmap.to_healpix(mask_iv,nside=aobj[qid].nside)
    #    w1[qid] = get_wfactor(mask_hp)[1]

    for rlz in tqdm.tqdm(aobj[qids[0]].rlz,desc='combine alms'):

        if misctools.check_path([aobj[qidc].falm['c']['T'][rlz],aobj[qidc].falm['n']['T'][rlz]],overwrite=overwrite,verbose=verbose): continue

        wTlm = {}
        for s in ['c','n']:
            for m in mtype:
                if rlz == 0 and s in ['n']: continue
                wTlm[s,m] = 0.
                for q in qids:
                    # walm = 1/N alm = 1/(N/r^2) alm/r
                    wTlm[s,m] += 1./ncl[q][:,None] * pickle.load(open(aobj[q].falm[s]['T'][rlz],"rb")) #/ w1[q]
                    #wTlm[s,m] += rcl[q][:,None]/(ncl[q][:,None]+ep) * pickle.load(open(aobj[q].falm[s]['T'][rlz],"rb")) #/ w1[q]
                    wTlm[s,m][:2,:] = 0.
        
                wTlm[s,m] *= Ncl[:,None]
                # save alm for each rlz
                pickle.dump((wTlm[s,m]),open(aobj[qidc].falm[s][m][rlz],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump((wTlm[s,m]),open(aobj[qidc].falm[s][m][rlz],"wb"),protocol=pickle.HIGHEST_PROTOCOL)


                

def diff_day_night(qid_d,qid_n,qid,overwrite=False,verbose=True,mtype=['T'],**kwargs):

    aobj = { q: local.init_analysis_params(qid=q,**kwargs) for q in [qid_d,qid_n,qid] }
    cl   = { s: np.zeros((len(aobj[qid].rlz),6,aobj[qid].lmax+1)) for s in ['c','s','n'] }

    for ii, rlz in enumerate(tqdm.tqdm(aobj[qid_d].rlz,desc='day - night')):

        if misctools.check_path(aobj[qid].falm['c']['T'][rlz],overwrite=overwrite,verbose=verbose): continue

        for s in ['c','s','n']:
            
            if rlz==0 and s in ['s','n']: continue

            # empirical Window correction from each day and night masks
            if s in ['s']:
                calm0 = pickle.load(open(aobj[qid_d].falm['c']['T'][rlz],"rb"))
                calm1 = pickle.load(open(aobj[qid_n].falm['c']['T'][rlz],"rb"))
                nalm0 = pickle.load(open(aobj[qid_d].falm['n']['T'][rlz],"rb"))
                nalm1 = pickle.load(open(aobj[qid_n].falm['n']['T'][rlz],"rb"))
                Talm0 = calm0 - nalm0
                Talm1 = calm1 - nalm1
            else:
                Talm0 = pickle.load(open(aobj[qid_d].falm[s]['T'][rlz],"rb"))
                Talm1 = pickle.load(open(aobj[qid_n].falm[s]['T'][rlz],"rb"))
            
            dTalm = Talm0 - Talm1

            if s == 'c':
                pickle.dump( dTalm, open(aobj[qid].falm[s]['T'][rlz],"wb"), protocol=pickle.HIGHEST_PROTOCOL )

            # aps
            cl[s][ii,0,:] = curvedsky.utils.alm2cl(aobj[qid].lmax,dTalm)
        
            # save cl for each rlz
            np.savetxt(aobj[qid].fcls[s][rlz],np.concatenate((aobj[qid].l[None,:],cl[s][ii,:,:])).T)

    # save mean cl to files
    if aobj[qid].rlz[-1] >= 2:
        if verbose:  print('save averaged diff day-night spectrum')
        imin = max(0,1-aobj[qid].rlz[0]) 
        for s in ['c','s','n']:
            np.savetxt(aobj[qid].fscl[s],np.concatenate((aobj[qid].l[None,:],np.mean(cl[s][imin:,:,:],axis=0),np.std(cl[s][imin:,:,:],axis=0))).T)

        

'''
def alm2aps_null(qid,overwrite=False,verbose=True,mtype=['T'],ep=1e-30,**kwargs):

    pid = qid.replace('_d0','_0')
    if 's16' in qid:
        pid = qid.replace('s16','boss')

    aobj = {q: local.init_analysis_params(qid=q,**kwargs) for q in [qid,pid]}

    mask_iv0 = load_mask(qid)
    mask_hp0 = enmap.to_healpix(mask_iv0,nside=aobj[qid].nside)
    mask_iv1 = load_mask(pid)
    mask_hp1 = enmap.to_healpix(mask_iv1,nside=aobj[pid].nside)

    w2 = np.average(mask_hp0**2)
    x2 = np.average(mask_hp0*mask_hp1)

    cl = np.zeros((len(aobj[qid].rlz),6,aobj[qid].lmax+1))
    xl = np.zeros((len(aobj[qid].rlz),6,aobj[qid].lmax+1))

    for ii, rlz in enumerate(tqdm.tqdm(aobj[qid].rlz)):

        #if misctools.check_path(aobj[qid].fnul['c'][rlz],overwrite=overwrite,verbose=verbose): continue

        Talm0 = pickle.load(open(aobj[qid].falm['c']['T'][rlz],"rb"))
        Talm1 = pickle.load(open(aobj[pid].falm['c']['T'][rlz],"rb"))
        cl[ii,0,:] = curvedsky.utils.alm2cl(aobj[qid].lmax,Talm0-Talm1) / w2
        xl[ii,0,:] = curvedsky.utils.alm2cl(aobj[qid].lmax,Talm0,Talm1) / x2
        
        # save cl for each rlz
        np.savetxt(aobj[qid].fcls_nul[rlz],np.concatenate((aobj[qid].l[None,:],cl[ii,:,:])).T)
        np.savetxt(aobj[qid].fcls_x[rlz],np.concatenate((aobj[qid].l[None,:],xl[ii,:,:])).T)

    # save mean cl to files
    if aobj[qid].rlz[-1] >= 2:
        if verbose:  print('cmb alm2aps: save sim')
        imin = max(0,1-aobj[qid].rlz[0]) 
        np.savetxt(aobj[qid].fscl_nul,np.concatenate((aobj[qid].l[None,:],np.mean(cl[imin:,:,:],axis=0),np.std(cl[imin:,:,:],axis=0))).T)
        np.savetxt(aobj[qid].fscl_x,np.concatenate((aobj[qid].l[None,:],np.mean(xl[imin:,:,:],axis=0),np.std(xl[imin:,:,:],axis=0))).T)
'''


