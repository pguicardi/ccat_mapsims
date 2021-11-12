import mapsims
import os
import sys
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
from pixell import enmap, enplot, reproject, utils, curvedsky
import pysm3.units as u
import pysm3 as pysm

def get_window(map,n_it=5):
    m_apo = map*0
    if hp.UNSEEN in map:
        m_apo[map!=hp.UNSEEN] = 1.0
    else:
        m_apo[map!=0] = 1.0
    i=0
    while i<= n_it:
        m_apo[m_apo<0.8] =0
        m_apo = hp.smoothing(m_apo, fwhm = 0.1)
        m_apo[m_apo<0] = 0
        m_apo /= np.max(m_apo)
        i+=1
    return m_apo

def lmax_f(NSIDE):
    return 3*NSIDE-1

def switch_to_local():
    (sys.path).insert(1,'/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/directories')
    return None

def transform_cmb_map(cmb, noise, ch):
    "ch is of the form 'LT6','LT5',etc... "
    chs = noise.tubes[ch]
    sky = pysm.Sky(nside=noise.nside, output_unit = cmb.input_units)
    sky.components.append(cmb)
    output_map = []
    for each in chs:
        print(each)
        bandpass_integrated_map = sky.get_emission(
            *each.bandpass
        ).value
        beam_width_arcmin = each.beam
        # smoothing and coordinate rotation with 1 spherical harmonics transform
        smoothed_map = hp.ma(
            pysm.apply_smoothing_and_coord_transform(
                bandpass_integrated_map,
                fwhm=beam_width_arcmin,
                lmax=3 * noise.nside - 1,
                rot=None,
                map_dist=None
            )
        )
        output_map.append(smoothed_map)
    output_map = np.array(output_map)
    return output_map

def apodize_map(map,n_itp =5):
    m_apo = get_window(map, n_it = n_itp)
    sh = m_apo.shape
    mt = m_apo.reshape(-1)
    npix = len(mt)
    non0 = np.where(mt!=0)
    w2 = np.sum(mt[non0]**2)/(npix)
    w4 = np.sum(mt[non0]**4)/(npix)
    mt = mt.reshape(sh)
    w = w2**2/w4
    tmp = map
    tmp1 = map==hp.UNSEEN
    tmp[tmp1] = 0
    return tmp*m_apo,w

def write_output_map(output,filename):
    for det in output.keys():
        hp.fitsfunc.write_map(filename[-5]+"_"+det+".fits", output[det], overwrite = True)
    return None

def read_output_map(filename,noise, ch):
    num = len(noise.tubes[ch[5:]])
    res = np.zeros((num,3,12*(noise.nside)**2))
    for i in np.arange(num):
        tmp = noise.tubes[ch[5:]][i]
        det = tmp.tag
        data = hp.fitsfunc.read_map(filename[-5]+"_"+det+".fits")
        res[i] = data

    return np.array(res)

def cls_to_dls(cls):
    ell = np.arange(np.array(cls).shape[-1])
    Dls = np.zeros(np.array(cls).shape)
    if len(cls.shape)>1:
        for i in np.arange(len(cls)):
            Dls[i] = cls[i]*ell*(ell+1) /2/np.pi
    else:
        Dls = cls*ell*(ell+1) /2/np.pi
    return Dls,ell

def calculate_power_spectrum(map,w, deconv_beam=False,beam_width = 1.0 ):
    cls = hp.sphtfunc.anafast(map)
    if deconv_beam:
        cls = deconvolve_beam(cls,beam_width)
    # Now we combine these to get the Dls
    Dls,ell = cls_to_dls(cls)
    return Dls/w,ell

def deconvolve_beam(cls, beam_width):
    """beam_width must be in degrees
       spec is the power spectrum to be deconvolved"""
            # beam_sig_rad = self.get_beams() * np.pi/180/60 / (8.*np.log(2))**0.5
            # beams = np.exp(-0.5 * ell*(ell+1) * beam_sig_rad[:,None]**2)
            # T_noise /= (beams[:,None,:] * beams[None,:,:])
            # P_noise /= (beams[:,None,:] * beams[None,:,:])
    ell = np.arange(cls.shape[-1])
    beam_sig_rad = beam_width*np.pi/180/60 / (8.*np.log(2))**0.5
    beam = np.exp(-0.5 * ell*(ell+1) * beam_sig_rad**2)
    result = np.array(cls)/(beam[None,:]*beam[None,:])
    return result[0]

def write_mapsims_sim(output_map):
    return None
