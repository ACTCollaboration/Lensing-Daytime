
# from external
import numpy as np
import healpy as hp
import sys
import pickle
import tqdm
import os    

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


def interface(run,qids,cqid=None,kwargs_ov={},kwargs={}):
    
    if 'alm' in run:
        if kwargs['fltr'] == 'cinv':
            sys.exit('interface function does not support cinv alm')
        else:
            map2alm(qids,**kwargs_ov,**kwargs)

    if 'aps' in run:
        alm2aps(qids,mtype=['T'],**kwargs_ov,**kwargs)

    if 'sup' in run:
        alm_supfac(qids,**kwargs_ov,**kwargs)

    if cqid is not None: # combine alms
        alm_comb(qids,cqid,mtype=['T'],**kwargs_ov,**kwargs)
        if kwargs['ivar'] == 'noivar': #otherwise need to speficy appropriate W factor for combined case
            alm2aps([cqid],mtype=['T'],**kwargs_ov,**kwargs)
        
    if 'each' in run:
        for qid_d, qid_n in zip(local.boss_d,local.boss_n):
            diff_day_night(qid_d,qid_n,'diff_'+qid_n,mtype=['T'],**kwargs_ov,**kwargs)

    if 's16' in run: # comb S16 day 
        print('combine alm for day')
        alm_comb(local.s_16_d,'comb_s16d',mtype=['T'],**kwargs_ov,**kwargs)    
        alm2aps(['comb_s16d'],mtype=['T'],**kwargs_ov,**kwargs)
        alm_supfac(['comb_s16d'],**kwargs_ov,**kwargs)

    if 'diff' in run:
    
        qid_dn = ['comb_d','comb_n']
        if kwargs['wind'] == 'iso15':
            if kwargs['ivar'] == 'V3':
                qids_d = ['boss_d02','boss_d03','boss_d04']
                qids_n = ['boss_02','boss_03','boss_04']
                #qids_d = local.boss_d
                #qids_n = local.boss_n
            else:
                qids_d = ['boss_d03']
                qids_n = ['boss_03','boss_04']
        elif kwargs['wind'] == 'com15':
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
        mask = load_survey_mask(qid,dg=1)
        
        # Define an object for sim generation
        model, season, array, patch, freq = local.qid_info(qid)
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

                # take coadd data
                maps['s'] = coadd_real_data(qid_array,mask) 

                # SZ map subtraction
                shape, wcs = maps['s'][0,0,:,:].shape, maps['s'].wcs
                if q in ['boss_d03','s16_d02','boss_03']: # need both freqs for combined data array of two freqs
                    sz_map_090 = enmap.project(enmap.read_map('data_local/input/S18d_202006_confirmed_model_f090.fits'),shape,wcs)/local.Tcmb
                    sz_map_150 = enmap.project(enmap.read_map('data_local/input/S18d_202006_confirmed_model_f150.fits'),shape,wcs)/local.Tcmb
                    maps['s'][0,0,:,:] -= sz_map_090
                    maps['s'][1,0,:,:] -= sz_map_150
                else:
                    sz_map = enmap.project(enmap.read_map('data_local/input/S18d_202006_confirmed_model_f150.fits'),shape,wcs)/local.Tcmb
                    maps['s'][0,0,:,:] -= sz_map                

            # simulation
            else:
                maps['s'], maps['n'], ivars = simobj.get_sim(season, patch, array, sim_num=i, set_idx=0)

                # coadd with civars weight implicitly multiplied
                civars = np.average(ivars,axis=1)
                civars[civars==0] = np.inf
                maps['s'] = mask[None,None,:,:] * np.average( maps['s']*ivars, axis=1 ) / civars / local.Tcmb 
                maps['n'] = mask[None,None,:,:] * np.average( maps['n']*ivars, axis=1 ) / civars / local.Tcmb

                # downgrade
                maps['s'] = enmap.downgrade(maps['s'], dg)
                maps['n'] = enmap.downgrade(maps['n'], dg)

            # save signal and noise to files
            for s in ['s','n']:
                if s=='n' and i==0: continue  # real data is saved to a signal file
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
    
    if aobj.ivar == 'vc':
        if '16' in aobj.wind:
            return hp.fitsfunc.read_map(aobj.fivar16,verbose=False)
        elif '15' in aobj.wind:
            return hp.fitsfunc.read_map(aobj.fivar15,verbose=False)
        else:
            sys.exit('ivar is not well defined')
    elif aobj.ivar == 'vd':
        ivar = enmap.read_map( aobj.fivarvd )
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    elif aobj.ivar == 'v3':
        ivar = enmap.read_map( local.init_analysis_params(qid='boss_d03').fivar )
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    elif aobj.ivar == 'V3':
        ivar = enmap.read_map( local.init_analysis_params(qid='boss_d03').fivar ) ** 2
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    elif aobj.ivar == 'noivar':
        return 1.
    elif aobj.ivar == 'base':
        ivar = enmap.read_map( aobj.fivar )
        return enmap.to_healpix(ivar[0],nside=aobj.nside)
    else:
        sys.exit('ivar is not well defined')
        

def create_ptsr_mask(nside,ascale=0.5,threshold=.0,radius=10.,ptsr='base'):

    # from catalogue
    #ras, decs, size = np.loadtxt('data_local/input/cat_crossmatched.txt',unpack=True,usecols=(0,1,5))
    #RAs  = ras[size>threshold]
    #DECs = decs[size>threshold]
    #arcm = size[size>threshold] * extend
    
    #RAs, DECs = interfaces.get_act_mr3f_union_sources(version='20200503_sncut_40')
    RAs, DECs = interfaces.get_act_mr3f_union_sources(version='20210209_sncut_10_aggressive') 
    arcm = np.ones(len(RAs)) * radius

    # additional mask
    if ptsr == 'PT':
        RAs_add  = np.array([187.3])
        DECs_add = np.array([2.])
        arcm_add = np.array([40.])
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
    # load combined window function in curvedsky (survey mask x inverse variance x point source mask)

    if survey_mask:
        # Load the default survey mask from soapack and project it onto healpix grids
        # (note that the mask is already apodized)

        if aobj.fltr=='cinv' or (aobj.fltr=='none' and aobj.wind=='base' and aobj.ascale==1.):
            # cinv case only loads survey_mask to make a boundary mask later
            if aobj.qid in local.qid_all: # for individual array, freqs
                mask = load_survey_mask(aobj.qid)
                mask = enmap.to_healpix(mask,nside=aobj.nside)  # to healpix map
            else:
                sys.exit('there are only survey masks for individual array maps')

        # Load other types of mask, including additional apodization (not equals to 1deg) and restricted boundary
        else:
            print('load custum mask',aobj.amask)
            mask = hp.fitsfunc.read_map( aobj.amask, verbose=False )
    else:
        mask = 1.

    # Load the default ivar map
    if with_ivar and aobj.ivar!='noivar' and aobj.fltr!='cinv':
        ivar = load_ivar_curvedsky(aobj)
        if ivar_norm:
            ivar = ivar/np.max(ivar)
    else:
        ivar = 1.

    # load ptsr mask for day
    if aobj.ptsr == '':
        ptsr = 1.
    else:
        print('load ptsr mask',aobj.fptsr)
        ptsr = hp.fitsfunc.read_map( aobj.fptsr ,verbose=False )

    return mask*ivar*ptsr


def get_wfactor(mask,wnmax=5):
    wn = np.zeros(wnmax)
    for n in range(1,wnmax):
        wn[n] = np.average(mask**n)
    Mask = mask*1.
    Mask[Mask!=0.] = 1.
    wn[0] = np.average(Mask)
    print('wfactors:',wn)
    return wn
    

def get_wfactors(qids,ascale,wind='base',ivar='base',ptsr='base',fltr='none',with_ivar=True,ivar_norm=False,wnmax=5):
    w = {}
    for q in qids:
        aobj = local.init_analysis_params(qid=q,ascale=ascale,wind=wind,ivar=ivar,ptsr=ptsr,fltr=fltr)
        mask = load_window_curvedsky(aobj,with_ivar=with_ivar,ivar_norm=ivar_norm)
        w[q] = get_wfactor(mask,wnmax=wnmax)
    return w


def beam_func(lmax,qid):
    
    dm = interfaces.models['dr5']()
    L  = np.arange(0, lmax+100, dtype=np.float)
    Bl = dm.get_beam_func(L, qid)
    return Bl[:lmax+1]

   
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


def remove_lxly(fmap,lmin=100,lmax=4096,lxcut=90,lycut=50):
    
    alm = enmap.fft(fmap)
    shape, wcs = fmap.shape, fmap.wcs
    kmask = define_lmask(shape,wcs,lmin=lmin,lmax=lmax,lxcut=lxcut,lycut=lycut)
    #print(np.shape(kmask),np.shape(alm))
    alm[kmask<0.5] = 0
    fmap_fl = enmap.ifft(alm).real
    return fmap_fl


def map2alm_core(fmap,aobj,amask):

    fmap_fl = remove_lxly( fmap, lmin=aobj.clmin, lmax=aobj.lmax )
    hpmap = enmap.to_healpix( fmap_fl, nside=aobj.nside )
    alm = curvedsky.utils.hp_map2alm( aobj.nside, aobj.lmax, aobj.lmax, hpmap*amask )

    return alm


def map2alm_core_spin(Qmap,Umap,aobj,amask):

    Qmap_fl = remove_lxly( Qmap, lmin=aobj.clmin, lmax=aobj.lmax )
    Umap_fl = remove_lxly( Umap, lmin=aobj.clmin, lmax=aobj.lmax )
    Qmap_hp = enmap.to_healpix( Qmap_fl, nside=aobj.nside )
    Umap_hp = enmap.to_healpix( Umap_fl, nside=aobj.nside )
    Ealm, Balm = curvedsky.utils.hp_map2alm_spin( aobj.nside, aobj.lmax, aobj.lmax, 2, Qmap_hp*amask, Umap_hp*amask )

    return Ealm, Balm


def map2alm(qids,overwrite=False,verbose=True,ascale=1.,**kwargs):

    # loop over qid
    for q in qids:
        
        aobj = local.init_analysis_params(qid=q,ascale=ascale,**kwargs)
        
        # beam
        Bl = beam_func(aobj.lmax,q)

        # total mask
        mask = load_window_curvedsky( aobj )

        # loop for map -> alm
        for i in tqdm.tqdm(aobj.rlz,desc='map -> alm'):

            for tp in ['T','P']:
                
                if tp=='T' and misctools.check_path(aobj.falm['c']['T'][i],overwrite=overwrite,verbose=verbose): continue
                if tp=='P' and misctools.check_path([aobj.falm['c']['E'][i],aobj.falm['c']['B'][i]],overwrite=overwrite,verbose=verbose): continue
            
                if i == 0: 

                    if tp=='T':
                        Tmap_c = enmap.read_map(aobj.fmap['s'][i])[0]
                        Talm_c = map2alm_core( Tmap_c, aobj, mask ) / Bl[:,None]
                        pickle.dump((Talm_c),open(aobj.falm['c']['T'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

                    if tp=='P':
                        Qmap_c = enmap.read_map(aobj.fmap['s'][i])[1]
                        Umap_c = enmap.read_map(aobj.fmap['s'][i])[2]
                        Ealm_c, Balm_c = map2alm_core_spin( Qmap_c, Umap_c, aobj, mask ) / Bl[:,None]
                        pickle.dump((Ealm_c),open(aobj.falm['c']['E'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                        pickle.dump((Balm_c),open(aobj.falm['c']['E'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

                else:

                    if tp=='T':
                        Tmap_s = enmap.read_map(aobj.fmap['s'][i])[0]
                        Tmap_n = enmap.read_map(aobj.fmap['n'][i])[0]
                        Talm_s = map2alm_core( Tmap_s, aobj, mask ) / Bl[:,None]
                        Talm_n = map2alm_core( Tmap_n, aobj, mask ) / Bl[:,None]
                        pickle.dump((Talm_n),open(aobj.falm['n']['T'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                        pickle.dump((Talm_s+Talm_n),open(aobj.falm['c']['T'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                        
                    if tp=='P':

                        Qmap_s = enmap.read_map(aobj.fmap['s'][i])[1]
                        Qmap_n = enmap.read_map(aobj.fmap['n'][i])[1]
                        Umap_s = enmap.read_map(aobj.fmap['s'][i])[2]
                        Umap_n = enmap.read_map(aobj.fmap['n'][i])[2]

                        Ealm_s, Balm_s = map2alm_core_spin( Qmap_s, Umap_s, aobj, mask ) / Bl[:,None]
                        Ealm_n, Balm_n = map2alm_core_spin( Qmap_n, Umap_n, aobj, mask ) / Bl[:,None]

                        pickle.dump((Ealm_n),open(aobj.falm['n']['E'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                        pickle.dump((Ealm_s+Ealm_n),open(aobj.falm['c']['E'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                        pickle.dump((Balm_n),open(aobj.falm['n']['B'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)
                        pickle.dump((Balm_s+Balm_n),open(aobj.falm['c']['B'][i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)



def alm2aps_core(lmax,alms,w2=1.):

    cl = np.zeros((6,lmax+1))
    
    mtype = alms.keys()
    
    # compute cls
    if 'T' in mtype:  
        #Talm = pickle.load(open(falm['T'],"rb"))
        cl[0,:] = curvedsky.utils.alm2cl(lmax,alms['T'])

    if 'E' in mtype:  
        #Ealm = pickle.load(open(falm['E'],"rb"))
        cl[1,:] = curvedsky.utils.alm2cl(lmax,alms['E'])

    if 'B' in mtype:  
        #Balm = pickle.load(open(falm['B'],"rb"))
        cl[2,:] = curvedsky.utils.alm2cl(lmax,alms['B'])

    if 'T' in mtype and 'E' in mtype:
        cl[3,:] = curvedsky.utils.alm2cl(lmax,alms['T'],alms['E'])

    if 'T' in mtype and 'B' in mtype:
        cl[4,:] = curvedsky.utils.alm2cl(lmax,alms['T'],alms['B'])

    if 'E' in mtype and 'B' in mtype:
        cl[5,:] = curvedsky.utils.alm2cl(lmax,alms['E'],alms['B'])

    return cl/w2
        

def alm2aps(qids,overwrite=False,verbose=True,mtype=['T','E','B'],cns=['c','n','s'],W2=None,**kwargs):

    for qid in qids: 
        
        aobj = local.init_analysis_params(qid=qid,**kwargs)
        if aobj.fltr == 'cinv': cns = ['c']
        
        if W2 is None:
            mask = load_window_curvedsky(aobj)
            w2 = get_wfactor(mask)[2]
        else:
            if isinstance(W2,dict):
                w2 = W2[qid]
            else:
                w2 = W2

        cl = {s: np.zeros((len(aobj.rlz),6,aobj.lmax+1)) for s in cns}

        for ii, rlz in enumerate(tqdm.tqdm(aobj.rlz,desc='alm -> aps')):

            if misctools.check_path(aobj.fcls['c'][rlz],overwrite=overwrite,verbose=verbose): continue

            if rlz == 0:
                alms = { m: pickle.load(open(aobj.falm['c'][m][rlz],"rb")) for m in mtype }
                cl['c'][ii,:,:] = alm2aps_core(aobj.lmax,alms,w2=w2)
            else:
                for s in cns:
                    if s=='s': continue # signal part will be computed from combined - noise
                    alms = { m: pickle.load(open(aobj.falm[s][m][rlz],"rb")) for m in mtype }
                    cl[s][ii,:,:]   = alm2aps_core(aobj.lmax,alms,w2=w2)

                # signal part
                if 's' in cns:
                    # combined - noise = signal
                    alms = { m: pickle.load(open(aobj.falm['c'][m][rlz],"rb")) - pickle.load(open(aobj.falm['n'][m][rlz],"rb")) for m in mtype }
                    cl['s'][ii,:,:] = alm2aps_cor(aobj.lmax,alms,w2=w2)
 
            # save cl for each rlz
            np.savetxt(aobj.fcls['c'][rlz],np.concatenate((aobj.l[None,:],cl['c'][ii,:,:])).T)
            if rlz >= 1:
                if 's' in cns:  np.savetxt(aobj.fcls['s'][rlz],np.concatenate((aobj.l[None,:],cl['s'][ii,:,:])).T)
                if 'n' in cns:  np.savetxt(aobj.fcls['n'][rlz],np.concatenate((aobj.l[None,:],cl['n'][ii,:,:])).T)

        # save mean cl to files
        if aobj.rlz[-1] >= 2:
            if verbose:  print('save averaged spectrum over rlz')
            imin = max(0,1-aobj.rlz[0]) 
            for s in cns:
                if misctools.check_path(aobj.fscl[s],overwrite=overwrite,verbose=verbose): continue
                np.savetxt(aobj.fscl[s],np.concatenate((aobj.l[None,:],np.mean(cl[s][imin:,:,:],axis=0),np.std(cl[s][imin:,:,:],axis=0))).T)


def comb_Nl(qids,ncl,rcl=None):

    Nl = 0.*ncl[qids[0]]
    
    for q in qids:
        if rcl is None:
            Nl[2:] += 1./ncl[q][2:]
        else:
            Nl[2:] += rcl[q][2:]**2/ncl[q][2:]
    Nl[2:] = 1./Nl[2:]
    return Nl



def alm_comb(qids,qidc,overwrite=False,verbose=True,mtypes=['T','E','B'],ep=1e-30,**kwargs):

    # qids = qid to be combined
    # qidc = output qid
    
    aobj = {q: local.init_analysis_params(qid=q,**kwargs) for q in qids+[qidc]}
    mid  = {'T':1,'E':2,'B':3}

    ncl = {}
    Ncl = {}
    for mi, m in mtype:
        # pre-computed noise spectra for inverse-variance weighting
        ncl[m]  = {q: (np.loadtxt(aobj[q].fscl['n'])).T[mid[m]] for q in qids}
        # norm
        Ncl[m]  = comb_Nl(qids,ncl)

    for rlz in tqdm.tqdm(aobj[qids[0]].rlz,desc='combine alms'):

        for m in mtype:
        
            if misctools.check_path([aobj[qidc].falm['c'][m][rlz],aobj[qidc].falm['n'][m][rlz]],overwrite=overwrite,verbose=verbose): continue

            walm = {}
            for s in ['c','n']:
            
                if rlz == 0 and s in ['n']: continue
                walm[s,m] = 0.
                for q in qids:
                    # walm = 1/N alm = 1/(N/r^2) alm/r
                    walm[s,m] += 1./ncl[m][q][:,None] * pickle.load(open(aobj[q].falm[s][m][rlz],"rb")) #/ w1[q]
                    #wTlm[s,m] += rcl[q][:,None]/(ncl[q][:,None]+ep) * pickle.load(open(aobj[q].falm[s]['T'][rlz],"rb")) #/ w1[q]
                    walm[s,m][:2,:] = 0.
        
                walm[s,m] *= Ncl[m][:,None]
                # save alm for each rlz
                pickle.dump((walm[s,m]),open(aobj[qidc].falm[s][m][rlz],"wb"),protocol=pickle.HIGHEST_PROTOCOL)


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



def wiener_cinv_core(qids,wqid,white=False,kwargs_ov={'overwrite':False,'verbose':True},kwargs_cmb={},kwargs_cinv={}):
    
    # parameter
    vb = kwargs_ov['verbose']
    
    # specify output filename
    aobj = local.init_analysis_params(qid=wqid,**kwargs_cmb)
    
    # input CMB data
    Aobj = {q: local.init_analysis_params(qid=q,**kwargs_cmb) for q in qids}
    
    mn  = len(qids)
    bl  = np.zeros((mn,aobj.lmax+1))
    Nij = np.zeros((1,mn,aobj.npix))
    mask = {}
    inl = np.ones((1,mn,aobj.lmax+1))
    
    for i, q in enumerate(qids):
    
        # beam
        bl[i,:] = beam_func(Aobj[q].lmax,Aobj[q].qid)

        # binary mask
        mask[q] = load_window_curvedsky( Aobj[q], with_ivar=False )
        mask[q][mask[q]!=0] = 1.
        
        # inv noise covariance
        #os.system('echo "'+q+',constructing noise covariance" >> test.txt')
        if vb:  print(q,'constructing noise covariance')

        Nvar = load_ivar_curvedsky(Aobj[q])

        if white:
            
            Nij[0,i,:] = mask[q] * Nvar/np.max(Nvar) / local.qid_wnoise(q)**2
            
        else:

            if q in local.boss_d:
                # averaged noise spectrum at S16 region
                bobj = local.init_analysis_params(qid=q,fltr='none',wind='com16',ivar='base',ptsr=aobj.ptsr)
            if q in local.boss_n or q in local.s_16_d:
                bobj = local.init_analysis_params(qid=q,fltr='none',wind='base',ivar='base',ptsr=aobj.ptsr)

            Nl  = np.loadtxt(bobj.fscl['n'],unpack=True)[1]
            inl[0,i,2:] = 1. / (Nl[2:]*bl[i,2:]**2)
            Nvar[Nvar<=0] = 1e-60
            Nij[0,i,:] = mask[q] * np.sqrt( Nvar/np.max(Nvar) ) #/ local.qid_wnoise(q)

        del Nvar

    # temperature map
    T  = np.zeros((1,mn,aobj.npix))
    
    for rlz in aobj.rlz:
    
        if misctools.check_path(aobj.falm['c']['T'][rlz],**kwargs_ov): continue
            
        #os.system('echo "'+str(rlz)+',loading temperature obs map" >> test.txt')
        if rlz==0: 
            for i, q in enumerate(qids):
                if vb:  print(rlz,q,'loading temperature obs map')
                Tc = enmap.read_map(Aobj[q].fmap['s'][rlz])[0]
                T[0,i,:] = mask[q] * enmap.to_healpix( remove_lxly( Tc, lmin=aobj.clmin, lmax=aobj.lmax ), nside=aobj.nside )
        else:
            for i, q in enumerate(qids):
                if vb:  print(rlz,q,'loading temperature sim map')
                Ts = enmap.read_map(Aobj[q].fmap['s'][rlz])[0]
                Tn = enmap.read_map(Aobj[q].fmap['n'][rlz])[0]
                T[0,i,:] = mask[q] * enmap.to_healpix( remove_lxly( Ts+Tn, lmin=aobj.clmin, lmax=aobj.lmax ), nside=aobj.nside )
    
        # cinv
        if vb:  print('cinv filtering')
        #os.system('echo "cinv filter" >> test.txt')
        if white:
            Tlm = curvedsky.cninv.cnfilter_freq(1,mn,aobj.nside,aobj.lmax,aobj.lcl[0:1,:],bl,Nij,T,verbose=kwargs_ov['verbose'],**kwargs_cinv)
        else:
            Tlm = curvedsky.cninv.cnfilter_freq(1,mn,aobj.nside,aobj.lmax,aobj.lcl[0:1,:],bl,Nij,T,verbose=kwargs_ov['verbose'],inl=inl,**kwargs_cinv)
        pickle.dump((Tlm[0,:,:]),open(aobj.falm['c']['T'][rlz],"wb"),protocol=pickle.HIGHEST_PROTOCOL)




# No longer used

'''

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
        rl = np.zeros((len(aobj.rlz),2,lmax+1))

        for ii, rlz in enumerate(tqdm.tqdm(aobj.rlz,desc='compute supfac')):
 
            # projected signal alm
            if aobj.fltr=='cinv':
                Tlm_c = pickle.load(open(aobj.falm['c']['T'][rlz],"rb"))
                Tlm_s = 0.*Tlm_c
            else:
                Tlm_c = pickle.load(open(aobj.falm['c']['T'][rlz],"rb"))
                Tlm_n = pickle.load(open(aobj.falm['n']['T'][rlz],"rb"))
                Tlm_s = Tlm_c - Tlm_n
            
            # input fullsky signal alm
            f = '/project/projectdirs/act/data/actsims_data/signal_v0.4/fullskyLensedUnabberatedCMB_alm_set00_'+str(rlz).zfill(5)+'.fits'
            Tlm_inp = np.complex128( hp.fitsfunc.read_alm( f, hdu = (1) ) ) / local.Tcmb
            ilmax = hp.sphtfunc.Alm.getlmax(len(Tlm_inp))
            Tlm_inp = curvedsky.utils.lm_healpy2healpix(Tlm_inp, ilmax)[:lmax+1,:lmax+1]

            # obs x input
            xl = curvedsky.utils.alm2cl(lmax,Tlm_c,Tlm_inp)
            sl = curvedsky.utils.alm2cl(lmax,Tlm_s,Tlm_inp)
            
            # input auto
            cl = curvedsky.utils.alm2cl(lmax,Tlm_inp)
            
            # take ratio to get suppression factor
            rl[ii,0,2:] = xl[2:]/cl[2:]
            rl[ii,1,2:] = sl[2:]/cl[2:]

        # save to file
        mrl = np.mean(rl,axis=0)
        np.savetxt( aobj.fsup, np.concatenate( ( aobj.l[None,:], mrl, mrl/w1 ) ).T )
        

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


def ivars_normalize(ivars):
    n1, n2, n3 = np.shape(ivars[:,:,:,0,0])
    nivars = 0.*ivars
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                nivars[i,j,k,:,:] = ivars[i,j,k,:,:]/np.max(ivars[i,j,k,:,:])
    return nivars


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


