
# from external
import numpy as np
import healpy as hp
import sys
import configparser
import pickle
from matplotlib.pyplot import *

# from act library
from pixell import enmap
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/soapack/")
from soapack import interfaces

# from cmblensplus/wrap/

# from cmblensplus/utils/
import misctools
import plottools as pl
import binning as bn
import analysis as ana


# fixed values
Tcmb = 2.726e6
ac2rad = np.pi/10800.

boss_d = ['boss_d01','boss_d02','boss_d03','boss_d04']
s_16_d = ['s16_d01','s16_d02','s16_d03']
boss_n = ['boss_01','boss_02','boss_03','boss_04']
boss_dn = boss_d + boss_n
day_all = boss_d + s_16_d
qid_all = day_all + boss_n

# qid for combined case
wqids = ['boss_s15d','boss_s15dn','boss_s15n','boss_s16d','boss_alld','boss_alldn']

#-----------------
# qid info
#-----------------

def get_subqids(wqid):
    
    if wqid == 'boss_s15d':  qids = boss_d
    if wqid == 'boss_s15n':  qids = boss_n
    if wqid == 'boss_s15dn': qids = boss_dn
    if wqid == 'boss_s16d':  qids = s_16_d
    if wqid == 'boss_alld':  qids = day_all
    if wqid == 'boss_alldn': qids = qid_all
    
    return qids


def qid_info(qid):

    if qid in ['bndn_01','bndn_02','boss_d01','boss_d02','boss_d03','boss_d04','s16_d01','s16_d02','s16_d03']:
        model = 'dr5'
    if qid in ['boss_01','boss_02','boss_03','boss_04','s16_01','s16_02','s16_03']:
        model = 'act_mr3'
    
    if model == 'act_mr3':

        if qid == 'boss_01':
            season, array, patch, freq = ('s15', 'pa1', 'boss','f150')
        if qid == 'boss_02':
            season, array, patch, freq = ('s15', 'pa2', 'boss','f150')
        if qid == 'boss_03':
            season, array, patch, freq = ('s15', 'pa3', 'boss','f090')
        if qid == 'boss_04':
            season, array, patch, freq = ('s15', 'pa3', 'boss','f150')

    if model == 'dr5':
        
        dm = interfaces.models['dr5']()
        season, array, patch, freq = (dm.ainfo(qid,'season'), dm.ainfo(qid,'array'), dm.ainfo(qid,'region'), dm.ainfo(qid,'freq'))
        
    return  model, season, array, patch, freq
     
    
def qid_label(qid):
    
    table = {'boss_01':'S15 PA1 f150 night', \
             'boss_02':'S15 PA2 f150 night', \
             'boss_03':'S15 PA3 f090 night', \
             'boss_04':'S15 PA3 f150 night', \
             'boss_d01':'S15 PA1 f150 day', \
             'boss_d02':'S15 PA2 f150 day', \
             'boss_d03':'S15 PA3 f090 day', \
             'boss_d04':'S15 PA3 f150 day', \
             's16_d01':'S16 PA2 f150 day', \
             's16_d02':'S16 PA3 f090 day', \
             's16_d03':'S16 PA3 f150 day', \
            }
    
    return  table[qid]


def qid_wnoise(qid):
    # white noise level at season 16 region (deepest region)
    
    table = {'boss_01':70., \
             'boss_02':35., \
             'boss_03':30., \
             'boss_04':45., \
             'boss_d01':70., \
             'boss_d02':40., \
             'boss_d03':30., \
             'boss_d04':45., \
             's16_d01':30., \
             's16_d02':30., \
             's16_d03':40., \
            }
    
    return  table[qid]*ac2rad/Tcmb


#------------------
# file directories
#------------------

# Define directory
def data_directory(root='/global/homes/t/toshiyan/Work/Ongoing/act_lens/'):
    
    direct = {}

    direct['root'] = root
    direct['input'] = '/project/projectdirs/act/data/actsims_data/signal_v0.4/'

    direct['dr4']  = root + 'data_dr4/'
    direct['mask'] = root + 'data_masks/'
    direct['plot'] = root + 'data_plots/'
    direct['sync'] = root + 'data_synced/'
    direct['local'] = root + 'data_local/'
    
    direct['cmb'] = direct['local'] + 'cmb/'

    return direct


#-----------------
# analysis object
#-----------------

# Define analysis parameters
class analysis_setup():

    def __init__(self,snmin=0,snmax=100,qid='boss_d01',fltr='none',lmin=1,lmax=4096,clmin=100,olmin=1,olmax=2048,bn=30,nside=2048,wind='base',ivar='base',ptsr='PT',ascale=1.):

        #//// load config file ////#
        conf = misctools.load_config('CMB')

        # rlz
        self.snmin = conf.getint('snmin',snmin)
        self.snmax = conf.getint('snmax',snmax)
        self.rlz   = np.linspace(self.snmin,self.snmax,self.snmax-self.snmin+1,dtype=np.int)
        if self.snmin == 0:
            self.snum  = self.snmax - self.snmin
        else:
            self.snum  = self.snmax - self.snmin + 1

        # multipole range of observed CMB alms
        self.lmin   = conf.getint('lmin',lmin)
        self.lmax   = conf.getint('lmax',lmax)
        
        # filtering multipole below clmin in addition to lx, ly before map->alm
        self.clmin  = conf.getint('clmin',clmin)

        # multipoles of output CMB spectrum
        self.olmin  = conf.getint('olmin',olmin)
        self.olmax  = conf.getint('olmax',olmax)
        self.bn     = conf.getint('bn',bn)
        self.binspc = conf.get('binspc','')

        # cmb map
        self.qid    = conf.get('qid',qid)
        self.fltr   = conf.get('fltr',fltr)

        # fullsky map
        self.nside  = conf.getint('nside',nside) #Nside for fullsky cmb map
        self.npix   = 12*self.nside**2

        # window, ivar, ptsr
        self.wind   = conf.get('wind',wind)
        if self.fltr == 'cinv':
            self.ascale = 0.
        else:
            self.ascale = conf.getfloat('ascale',ascale)
        self.apotag = 'a'+str(self.ascale)+'deg'
        self.ivar   = conf.get('ivar',ivar)
        self.ptsr   = conf.get('ptsr',ptsr)

        if self.fltr == 'cinv':
            self.wtype = '_'.join( [ self.wind, self.ivar, self.ptsr ] )
        else:
            self.wtype = '_'.join( [ self.wind, self.ivar, self.ptsr, self.apotag ] )

        # do
        self.doreal = conf.getboolean('doreal',False)
        self.dodust = conf.getboolean('dodust',False)


    def filename(self):

        #//// root directories ////#
        d = data_directory()
        d_map = d['cmb'] + 'map/'
        d_alm = d['cmb'] + 'alm/'
        d_aps = d['cmb'] + 'aps/'
        d_msk = d['cmb'] + 'mask/'

        #//// index ////#
        self.ids = [str(i).zfill(5) for i in range(-1,1000)]
        self.ids[0] = 'real'  # change 1st index
        ids = self.ids

        #//// Partial sky CMB maps from actsims ////#
        self.fmap = { s: [d_map+s+'_'+self.qid+'_'+x+'.fits' for x in ids] for s in ['s','n'] }

        # ivar maps
        self.fivar   = d_map+'ivar_'+self.qid+'.fits'
        self.fivar15 = d_map+'ivar_com15.fits'
        self.fivar16 = d_map+'ivar_com16.fits'
        self.fivarvd = d_map+'ivar_comdy.fits'

        # input klm realizations
        self.fiklm = [ d['input'] + 'fullskyPhi_alm_'+x+'.fits' for x in ids[1:] ]

        #//// base best-fit cls ////#
        # aps of Planck 2015 best fit cosmology
        self.fucl = d['local'] + 'input/cosmo2017_10K_acc3_scalCls.dat'
        self.flcl = d['local'] + 'input/cosmo2017_10K_acc3_lensedCls.dat'

        #//// basic tags for alm ////#
        self.stag = '_'.join( [ self.qid , self.wtype , self.fltr , 'lc'+str(self.clmin) ] )
        self.ntag = '_'.join( [ self.qid , self.wtype , self.fltr , 'lc'+str(self.clmin) ] )

        #//// Derived data filenames ////#
        # cmb signal/noise alms
        self.falm, self.fscl, self.fcls = {}, {}, {}
        for s in ['s','n','p','c']:
            if s in ['n','p']: tag = self.ntag
            if s in ['s','c']: tag = self.stag
            self.falm[s] = { m: [d_alm+'/'+s+'_'+m+'_'+tag+'_'+x+'.pkl' for x in ids] for m in ['T','E','B'] }

            # cmb aps
            self.fscl[s] = d_aps+'aps_sim_1d_'+tag+'_'+s+'.dat'
            self.fcls[s] = [d_aps+'/rlz/cl_'+tag+'_'+s+'_'+x+'.dat' for x in ids]

        # null cmb aps
        self.fscl_nul = d_aps+'aps_sim_1d_'+self.ntag+'_null.dat'
        self.fcls_nul = [d_aps+'/rlz/cl_'+self.ntag+'_null_'+x+'.dat' for x in ids]
        self.fscl_x = d_aps+'aps_sim_1d_'+self.stag+'_cross.dat'
        self.fcls_x = [d_aps+'/rlz/cl_'+self.stag+'_cross_'+x+'.dat' for x in ids]
        
        # suppression
        self.fsup = d_aps + 'supfac_'+self.stag+'.dat'

        # beam
        self.fbeam = d['local'] + 'beam/' + self.qid+'.dat'
        
        # custom mask
        if self.wind == 'base': # no restriction to area
            self.amask = d_msk + self.qid+'_base_'+self.apotag+'.fits'
        else:
            self.amask = d_msk + self.wind+'_'+self.apotag+'.fits'
        
        # ptsr mask
        #self.fptsr_old = d_msk + 'custom_ptsr_square_mask.fits'
        #self.fptsr = d_msk + 'ptsr_cat_crossmatched.fits'
        self.fptsr = d_msk + 'ptsr_'+self.ptsr+'.fits'

        #//// basic tags ////#
        # output multipole range
        self.otag = '_oL'+str(self.olmin)+'-'+str(self.olmax)+'_b'+str(self.bn)


    def array(self):

        #multipole
        self.l  = np.linspace(0,self.lmax,self.lmax+1)
        self.kL = self.l*(self.l+1)*.5

        #theoretical cl
        self.ucl = np.zeros((5,self.lmax+1)) # TT, EE, TE, pp, Tp
        self.ucl[:,2:] = np.loadtxt(self.fucl,unpack=True,usecols=(1,2,3,4,5))[:,:self.lmax-1] 
        self.ucl[:3,:] *= 2.*np.pi / (self.l**2+self.l+1e-30) / Tcmb**2
        self.ucl[3,:] *= 1. / (self.l+1e-30)**4 / Tcmb**2

        self.lcl = np.zeros((4,self.lmax+1)) # TT, EE, BB, TE
        self.lcl[:,2:] = np.loadtxt(self.flcl,unpack=True,usecols=(1,2,3,4))[:,:self.lmax-1] 
        self.lcl *= 2.*np.pi / (self.l**2+self.l+1e-30) / Tcmb**2
        
        self.cpp = self.ucl[3,:]
        self.ckk = self.ucl[3,:] * (self.l**2+self.l)**2/4.


#----------------
# initial setup
#----------------

def init_analysis_params(**kwargs):
    # setup parameters, filenames, and arrays
    aobj = analysis_setup(**kwargs)
    analysis_setup.filename(aobj)
    analysis_setup.array(aobj)
    return aobj
      

# ----------
# Test
# ----------

def quick_rec(alm,ocl,lcl,al,mask=None,rlmin=500,rlmax=3000,nside=2048,lmax=2048):
    import curvedsky
    if mask is None:
        wTlm = alm[:rlmax+1,:rlmax+1]
    else:
        wTlm = curvedsky.utils.mulwin(alm[:rlmax+1,:rlmax+1],mask)
    fTlm = wTlm/(ocl[:rlmax+1,None]+1e-30)
    klm, __ = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,lcl[:rlmax+1],fTlm,fTlm,gtype='k',nside_t=nside)
    klm *= al[:lmax+1,None]
    kl = curvedsky.utils.alm2cl(lmax,klm)
    return klm, kl


# ----------
# Plot
# ----------

def show_tmap(Tlm,ocl,mask=1,lmin=500,lmax=3000,v=3e11,nside=512,lonra=[148,243],latra=[-3,20],title=''):
    import curvedsky
    Flm = Tlm.copy()
    Flm[:lmin,:] = 0.
    Tmap = curvedsky.utils.hp_alm2map(nside,lmax,lmax,Flm[:lmax+1,:lmax+1]/(ocl[:lmax+1,None]+1e-30))
    hp.cartview(mask*Tmap,lonra=lonra,latra=latra,min=-v,max=v,cbar=False,title=title)

    
def show_kmap(klm=None,fname=None,lmin=200,lmax=1024,nside=1024,v=.1,lonra=[147,244],latra=[-3,21],output=False,title=''):
    import curvedsky
    if fname is not None:
        Flm, __ = pickle.load(open(fname,"rb"))
    if klm is not None:
        Flm = klm.copy()
    Flm[:lmin,:] = 0.
    kmap = curvedsky.utils.hp_alm2map(nside,lmax,lmax,Flm[:lmax+1,:lmax+1])
    hp.cartview(kmap,lonra=lonra,latra=latra,min=-v,max=v,cbar=False,title=title)
    if output:
        return kmap


def load_spec(qobj,mb,rlz=None,cn=1,outN0=False):
    # load data
    l, al = (np.loadtxt(qobj.f['TT'].al,usecols=(0,cn))).T
    l, n0 = (np.loadtxt(qobj.f['TT'].n0bs,usecols=(0,cn))).T
    if rlz is None:
        #rd = n0.copy()
        rd = (np.loadtxt(qobj.f['TT'].rdn0[0])).T[cn]
        fcl = qobj.f['TT'].cl[:101]
    else:
        rd = np.array( [ (np.loadtxt(qobj.f['TT'].rdn0[i])).T[cn] for i in rlz ] )
        fcl = [ qobj.f['TT'].cl[i] for i in rlz ]
    # binning 
    vl = al/np.sqrt(l+.5+1e-30)
    nb = bn.binning(n0,mb,vl=vl)
    rb = bn.binning(rd,mb,vl=vl)
    mkk, __, skk, okk = bn.binned_spec(mb,fcl,cn=cn,doreal=True,opt=True,vl=vl)
    # obs and sim kk
    if rlz is None:
        Okk = okk-rb-mkk/100.
        Skk = skk-nb-mkk/99.
    else:
        Okk = okk-rb[0]-mkk/100.
        Skk = skk-rb[1:,:]-mkk/99.
    # mean and var of sim kk
    Mkk = np.mean(Skk,axis=0)
    Vkk = np.std(Skk,axis=0)
    # output
    if outN0:
        return Mkk, Vkk, Skk, Okk, nb
    else:
        return Mkk, Vkk, Skk, Okk


def plot_spec_kk(qobj,rlz=None,cn=1,lmin=40,lmax=2048,bnum=10,output=True,verbose=True,lfac=0.0,plot_real=False,yrange=False):
    # compute binned spectrum
    mb = bn.multipole_binning(bnum,lmin=lmin,lmax=lmax)
    Mkk, Vkk, Skk, Okk, nb = load_spec(qobj,mb,rlz=rlz,cn=cn,outN0=True)
    if verbose: print(np.sqrt(np.sum(Mkk**2/Vkk**2)))
    # statistics
    st = ana.amplitude(Mkk,Skk,fcb=None,diag=False,disp=True)
    # plot
    if yrange:
        pl.plot_1dstyle(fsize=[10,4],xmin=mb.lmin,xmax=lmax,ymin=-2e-8,ymax=1.2e-7,ylab='$L^{'+str(lfac)+'}C_L^{\kappa\kappa}$')
    else:
        pl.plot_1dstyle(fsize=[10,4],xmin=mb.lmin,xmax=lmax,ylab='$L^{'+str(lfac)+'}C_L^{\kappa\kappa}$')
    aobj = init_analysis_params()
    s = mb.bc**lfac
    if plot_real:
        errorbar(mb.bc+5,s*Okk,yerr=s*Vkk,fmt='o')
    else:
        errorbar(mb.bc+5,s*Mkk,yerr=s*Vkk,fmt='o')
        plot(aobj.l,aobj.l**lfac*aobj.ckk,color='k',ls='--')
    axhline(0,color='k')
    show()
