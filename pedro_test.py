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

ch = "tube:LC1"

# noise = mapsims.SONoiseSimulator(
#     nside=NSIDE,
#     return_uK_CMB = True,
#     sensitivity_mode = "baseline",
#     apply_beam_correction = False,
#     apply_kludge_correction = True,
#     rolloff_ell = 50,
#     survey_efficiency = 0.2,
#     full_covariance = True,
#     LA_years = 5,
#     LA_noise_model = "SOLatV3point1",
#     elevation = 50,
#     SA_years = 5,
#     SA_one_over_f_mode = "pessimistic"
# )

noise = mapsims.SONoiseSimulator(
    nside=NSIDE,
    return_uK_CMB = True,
    sensitivity_mode = "baseline",
    apply_beam_correction = False,
    apply_kludge_correction = False,
    homogeneous=True,
    rolloff_ell = 50,
    survey_efficiency = 1.0,
    full_covariance = False,
    LA_years = 5,
    LA_noise_model = "CcatLatv2b",
    elevation = 50,
    SA_years = 5,
    SA_one_over_f_mode = "pessimistic"
)


pysm_string = "none"
tag = pysm_string

filename = ch[5:]+"_NSIDE_" + str(NSIDE) + "_TAG_" +"hom_"+ tag + ".fits"
simulator = mapsims.MapSim(
    channels=ch,
    nside=NSIDE,
    unit="uK_CMB",
    pysm_output_reference_frame="C",
#    pysm_components_string=pysm_string,
    output_filename_template = filename,
    pysm_custom_components={"cmb": cmb},
    other_components={"noise": noise},
)
output_map_full = simulator.execute()
write_output_map(output_map_full,"output/"+filename)

output_map_full = read_output_map("output/"+filename,noise, ch)
#noise_sim = noise.simulate(ch[5:], nsplits = 1)

for i in np.arange(2):

    output_T_det1= output_map_full[i][0]

    apodized_map,w = apodize_map(output_T_det1, n_itp=5)
    tmp = get_window(output_T_det1, n_it =5)

    hp.mollview(apodized_map)
    plt.show()
    plt.clf()

    #Compute PS for CMB and Total Map
    Dl_cmb, ell_c = cls_to_dls(hp.alm2cl(cmb.alm))
    Dl_apo,ell = calculate_power_spectrum(apodized_map,w,deconv_beam=True,beam_width = noise.get_beam_fwhm(ch[5:])[i])

    #Compute Noise sims

#    apo_noise1,w = apodize_map(noise_sim[i][0][0])
#    noise_spec1,elln = calculate_power_spectrum(apo_noise1,w)

    #Compute PS for Noise Channel
    ell_noise,cls_T,cls_P = noise.get_fullsky_noise_spectra(ch[5:])
    Dl_noise1 = (cls_T[i])*ell_noise*(ell_noise+1)/(2*np.pi)

    if True: #PLOT CMB SPECS
        plt.plot(ell_c, Dl_cmb[0], '.', label = 'Theory CMB')
        plt.plot(ell, Dl_apo,'.', label = 'Deconvolved Simulated CMB + Noise')
        plt.plot(ell_noise, Dl_noise1, '.', label = 'Theoretical Noise Curve')
        plt.legend()
        plt.ylabel(r'$D_\ell$')
        plt.xlabel(r'$\ell$')
        #plt.xscale("log")
        plt.show()
