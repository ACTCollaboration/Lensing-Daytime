
# from external
import numpy as np
import healpy as hp
import sys
import configparser
import pickle

# from act library
from pixell import enmap
sys.path.append("/global/homes/t/toshiyan/Work/Lib/actlib/soapack/")
from soapack import interfaces

# from cmblensplus/wrap/
import curvedsky
import basic

# from cmblensplus/utils/
import constants
import misctools


# fixed values
Tcmb = 2.726e6


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


# Define analysis parameters
class analysis_setup():

    def __init__(self,snmin=0,snmax=100,qid='boss_d01',fltr='none',lmin=1,lmax=4096,olmin=1,olmax=2048,bn=30,nside=2048,wtype='',ascale=1.):

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

        # window
        self.wtype  = conf.get('wtype',wtype)
        self.ascale = conf.getfloat('ascale',ascale)

        # do
        self.doreal = conf.getboolean('doreal',False)
        self.dodust = conf.getboolean('dodust',False)


    def filename(self):

        #//// root directories ////#
        d = data_directory()
        d_map = d['cmb'] + 'map/'
        d_alm = d['cmb'] + 'alm/'
        d_aps = d['cmb'] + 'aps/'

        #//// basic tags ////#
        # for alm
        apotag = 'a'+str(self.ascale)+'deg'
        self.stag = '_'.join( [ self.qid , self.wtype , apotag , self.fltr ] )
        self.ntag = '_'.join( [ self.qid , self.wtype , apotag , self.fltr ] )

        # output multipole range
        self.otag = '_oL'+str(self.olmin)+'-'+str(self.olmax)+'_b'+str(self.bn)

        #//// index ////#
        self.ids = [str(i).zfill(5) for i in range(-1,1000)]
        self.ids[0] = 'real'  # change 1st index
        ids = self.ids

        #//// Partial sky CMB maps from actsims ////#
        self.fmap = {}
        for s in ['s','n']:
            self.fmap[s] = [d_map+s+'_'+self.qid+'_'+x+'.fits' for x in ids]

        self.fivar = d_map+'ivar_'+self.qid+'.fits'

        # input klm realizations
        self.fiklm = [ d['input'] + 'fullskyPhi_alm_'+x+'_klm.fits' for x in ids ]

        #//// base best-fit cls ////#
        # aps of Planck 2015 best fit cosmology
        self.fucl = d['local'] + 'input/cosmo2017_10K_acc3_scalCls.dat'
        self.flcl = d['local'] + 'input/cosmo2017_10K_acc3_lensedCls.dat'

        #//// Derived data filenames ////#
        # cmb signal/noise alms
        self.falm = {}
        self.fscl = {}
        self.fcls = {}
        for s in ['s','n','p','c']:
            self.falm[s] = {}
            if s in ['n','p']: tag = self.ntag
            if s in ['s','c']: tag = self.stag
            for m in ['T','E','B']:
                self.falm[s][m] = [d_alm+'/'+s+'_'+m+'_'+tag+'_'+x+'.pkl' for x in ids]

            # cmb aps
            self.fscl[s] = d_aps+'aps_sim_1d_'+tag+'_'+s+'.dat'
            self.fcls[s] = [d_aps+'/rlz/cl_'+tag+'_'+s+'_'+x+'.dat' for x in ids]

        # null cmb aps
        self.fscl_nul = d_aps+'aps_sim_1d_'+tag+'_null.dat'
        self.fcls_nul = [d_aps+'/rlz/cl_'+tag+'_null_'+x+'.dat' for x in ids]
        self.fscl_x = d_aps+'aps_sim_1d_'+tag+'_cross.dat'
        self.fcls_x = [d_aps+'/rlz/cl_'+tag+'_cross_'+x+'.dat' for x in ids]

        # beam
        self.fbeam = d['local'] + 'beam/' + self.qid+'.dat'


    def array(self):

        #multipole
        self.l  = np.linspace(0,self.lmax,self.lmax+1)
        self.kL = self.l*(self.l+1)*.5

        #binned multipole
        self.bp, self.bc = basic.aps.binning(self.bn,[self.olmin,self.olmax],self.binspc)

        #theoretical cl
        self.ucl = np.zeros((5,self.lmax+1)) # TT, EE, BB, TE
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
      

def qid_info(qid):

    if qid in ['boss_d01','boss_d02','boss_d03','boss_d04','s16_d01','s16_d02','s16_d03']:
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
     
    


