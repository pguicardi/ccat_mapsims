import sys
(sys.path).insert(1,'/Users/pedroguicardi/Desktop/CMB_Analysis/MAPSIMS/directories')

import mapsims
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import numpy as np
import healpy as hp
import matplotlib
from matplotlib import pyplot as plt
from pixell import enmap, enplot, reproject, utils, curvedsky
from ad_fns import *
from astropy.io import fits
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})




NSIDE = 256
ell_cutoff = False
cmb = mapsims.SOPrecomputedCMB(
    num=2,
    nside=NSIDE,
    lensed=False,
    aberrated=False,
    has_polarization=True,
    cmb_set=0,
    cmb_dir="mapsims/mapsims/tests/data",
    input_units="uK_CMB",
)

noise = mapsims.SONoiseSimulator(
    nside=NSIDE,
    return_uK_CMB = True,
    sensitivity_mode = "baseline",
    apply_beam_correction = False,
    apply_kludge_correction = True,
    rolloff_ell = 50,
    survey_efficiency = 0.2,
    full_covariance = True,
    LA_years = 5,
    LA_noise_model = "CcatLatv2b",
    elevation = 50,
    SA_years = 5,
    SA_one_over_f_mode = "pessimistic"
)

noise2 = mapsims.SONoiseSimulator(
    nside=NSIDE,
    return_uK_CMB = True,
    sensitivity_mode = "baseline",
    apply_beam_correction = False,
    apply_kludge_correction = True,
    rolloff_ell = 50,
    survey_efficiency = 0.2,
    full_covariance = True,
    LA_years = 5,
    LA_noise_model = "SOLatV3point1",
    elevation = 50.0,
    SA_years = 5,
    SA_one_over_f_mode = "pessimistic"
)

ch = "tube:LC1"

#Compute PS for Noise Channel
if True: #Plot Noise Spectra
    ell_noise,cls_T,cls_P = noise.get_fullsky_noise_spectra(ch[5:])
    #Dl_noise1 = (cls_T[0])*ell_noise*(ell_noise+1)/(2*np.pi*w)
    #Dl_noise2 = (cls_T[1])*ell_noise*(ell_noise+1)/(2*np.pi*w)

    plt.plot(ell_noise, cls_T[0], label = str(noise.tubes[ch[5:]][0]) )
    plt.plot(ell_noise, cls_T[1], label = str(noise.tubes[ch[5:]][1]) )
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([10**2, 10**4])
    plt.ylim([10**(-5), 10**3])
    plt.ylabel(r'$C_{\ell}$')
    plt.xlabel(r'$\ell$')
    #plt.xscale("log")
    plt.show()



pysm_string = "d0"
tag = pysm_string

filename = ch[5:]+"_NSIDE_" + str(NSIDE) + "_TAG_" + tag + ".fits"

simulator = mapsims.MapSim(
    channels=ch,
    nside=NSIDE,
    unit="uK_CMB",
    pysm_output_reference_frame="C",
    pysm_components_string=pysm_string,
    output_filename_template = filename,
    pysm_custom_components={"cmb": cmb},
    other_components={"noise": noise},
)

output_map_full = simulator.execute()
write_output_map(output_map_full,"output/"+filename)

output_map_full = read_output_map("output/"+filename,noise, ch)
output_T_det1= output_map_full[1][0]
#output_T_det2 = output_map_full[1][0]
# plot_maps = False
# tmp = [output_T_det2]
# if plot_maps:
#     hp.mollview(tmp[0], max = 300, min = -300, unit="uK_CMB", title = "T noise map ")
#     #hp.mollview(tmp[1], max = 30, min = -30, unit="uK_CMB", title = "Q noise map ")
#     #hp.mollview(tmp[2], max = 30, min = -30, unit="uK_CMB", title = "U noise map ")
#     plt.show()
#     plt.clf()



# tmp = get_window(output_T_det1, n_it =5)
# hp.mollview(tmp, max = 1, min = -1, unit="uK_CMB", title = "Apodized T ")
# plt.show()
# plt.clf()

apodized_map,w = apodize_map(output_T_det1, n_it=5)
hp.mollview(output_T_det1, max = 300, min = -300, unit="uK_CMB", title = "Apodized T ")
plt.show()
plt.clf()

exit()

#Compute PS for CMB and Total Map
Dl_cmb, ell_c = cls_to_dls(hp.alm2cl(cmb.alm))
Dl,ell = calculate_power_spectrum(output_T_det1,w)
Dl_apo,ell = calculate_power_spectrum(apodized_map,w)

#Compute PS for Noise Channel
ell_noise,cls_T,cls_P = noise.get_fullsky_noise_spectra(ch[5:])
Dl_noise1 = (cls_T[0])*ell_noise*(ell_noise+1)/(2*np.pi)
Dl_noise2 = (cls_T[1])*ell_noise*(ell_noise+1)/(2*np.pi)


plt.plot(ell_noise, Dl_noise2, '.', label = 'Noise ')
# plt.plot(ell_noise, Dl_noise2, '.', label = 'Noise ')
plt.plot(ell_c, Dl_cmb[0], '.', label = 'CMB Map')
plt.plot(ell, Dl_apo, '.', label = 'Apodized Map ')
#plt.plot(ell, Dl, '.', label = 'Non-Apodized Map')
plt.legend()
plt.ylabel("D_ell")
plt.xlabel('ell')
#plt.xscale("log")
plt.show()
