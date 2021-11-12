import numpy as np
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import healpy as hp
import matplotlib
from matplotlib import pyplot as plt
from ad_fns import *
import healpy as hp


def get_Dls(filename,filepath,cols_to_read):
    """cols_to_read must be a tuple"""
    data = np.loadtxt(filepath+filename, usecols=cols_to_read, unpack=True)
    ell = data[0]
    dls = data[1:]

    return ell,dls

def compute_cls(ell,dls):
    cls = np.zeros(dls.shape)
    for i in np.arange(dls.shape[0]):
        cls[i] = dls[i]*2 * np.pi / (ell*(ell+1.))
        cls[i][0] = 0
        cls[i][1] = 0
    return np.array(cls)


NSIDE =256
n = 2 #number to store as

file_dir = "/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/mapsims/mapsims/tests/data"
filename_scal = "/camb_08284993_scalcls.fits"
cls = hp.fitsfunc.read_cl(file_dir+filename_scal)*10**12 #convert to uK
alms = hp.sphtfunc.synalm(cls,new=True)
Dls_T,ell = cls_to_dls(cls[0])



# ordering = ["TT","EE","BB", "TE"]
# cls=compute_cls(ell,Dls)

# alm = hp.synalm(tt)
# m = hp.alm2map(alm,nside)
#
# cl = hp.anafast(m)
# cl_ell = np.arange(len(cl))

if False:
    map = hp.alm2map(alms,NSIDE)
    hp.mollview(map[0], max = 300, min = -300, unit="uK_CMB", title = "CMB I map ")
    plt.show()
    plt.clf()
    hp.mollview(map[1],  unit="uK_CMB", title = "CMB Q map ")
    plt.show()
    plt.clf()
    hp.mollview(map[2],  unit="uK_CMB", title = "CMB U map ")
    plt.show()
    plt.clf()

try:
    hp.fitsfunc.write_alm("/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/mapsims/mapsims/tests/data/fullskyUnlensedUnabberatedCMB_alm_set00_0000" +str(n)+".fits", alms)
except:
    os.remove("/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/mapsims/mapsims/tests/data/fullskyUnlensedUnabberatedCMB_alm_set00_0000" +str(n)+".fits")
    hp.fitsfunc.write_alm("/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/mapsims/mapsims/tests/data/fullskyUnlensedUnabberatedCMB_alm_set00_0000" +str(n)+".fits", alms)

alms = hp.read_alm("/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/mapsims/mapsims/tests/data/fullskyUnlensedUnabberatedCMB_alm_set00_0000" +str(n)+".fits")
cls = hp.sphtfunc.alm2cl(alms)

Dls_a,ell_a = cls_to_dls(cls)
print(alms.shape, cls.shape, Dls_a.shape)
if True:
    plt.plot(ell, Dls_T, '.', label = 'Temperature B4 map')
    plt.plot(ell_a, Dls_a, '.', label = 'Temperature After')
    #plt.plot(ell, Dl, '.', label = 'Non-Apodized Map')
    plt.legend()
    plt.ylabel("D_ell")
    plt.xlabel('ell')
    #plt.xscale("log")
    plt.show()
