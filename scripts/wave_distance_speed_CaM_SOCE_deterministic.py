import os
import h5py
import numpy as np
import sys
import argparse
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
import utility_functions as utils



colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}

stim_dend = "dend26"

directories = {
    "Ca_wave_simple_SERCA_SOCE_det_lim": "model_%s%s_simple_SERCA_%s%s_diam_%s_um_50_um_%s_nM.h5",
     #"Ca_wave_simple_SERCA_SOCE_10xRyR_det_lim": "model_%s%s_simple_SERCA_%s%s_diam_%s_um_50_um_%s_nM.h5",
    "Ca_wave_simple_SERCA_SOCE": "model_%s%s_simple_SERCA_%s%s_diam_%s_um_50_um_%s_nM.h5",
    

    
}

descr = {
    "Ca_wave_simple_SERCA_SOCE": "RyR",
    "Ca_wave_simple_SERCA_SOCE_10xRyR_det_lim":"RyR",
    "Ca_wave_simple_SERCA_SOCE_det_lim": "RyR",

}   

types = {
    "Ca_wave_simple_SERCA_SOCE": "stochastic",
    "Ca_wave_simple_SERCA_SOCE_det_lim": "deterministic limit",
    "Ca_wave_simple_SERCA_SOCE_10xRyR_det_lim":"deterministic limit 2",
}

marker = {
    "Ca_wave_simple_SERCA_SOCE": "full",
    "Ca_wave_simple_SERCA_SOCE_det_lim": "none",
    "Ca_wave_simple_SERCA_SOCE_10xRyR_det_lim":"none",
}

stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]

symbol = {
    "Ca_wave_simple_SERCA_SOCE": "o",
    "Ca_wave_simple_SERCA_SOCE_det_lim": "o",
    "Ca_wave_simple_SERCA_SOCE_10xRyR_det_lim":"d",
 }


    
NA = Avogadro*1e-23
specie_dict = {
    "Ca": ["Ca"],
    "CaOut": ["CaOut"],
    "CaER": ["CaER"],
    "RyRO": ["RyRO1", "RyRO2", "RyRCaMO1", "RyRCaMO2"],
    "STIM_CaER": ["STIM_2CaER"],
    "Orai": ["OraiSTIM_4", "Orai2STIM_4", "Orai3STIM_4"]
}
multiplier = {
    "Ca": 1,
    "CaOut": 1,
    "CaER": 1,
    "RyRO1": 1,
    "RyRO2": 1,
    "RyRCaMO1": 1,
    "RyRCaMO2": 1,
    "STIM_2CaER": 1,
    "OraiSTIM_4": 1,
    "Orai2STIM_4": 2,
    "Orai3STIM_4": 3,
}
label_list = []
def Parser():
    parser = argparse.ArgumentParser(description='Generate figs of avg conc')
    parser.add_argument('--species', default="Ca",
                        help='Ca, RyRO, CaER, CaOut, RyR')

    return parser


s = ""


if __name__ == '__main__':
    fnames = []
    args = Parser().parse_args()
    chosen_specie = args.species
    if chosen_specie in ["Ca", "CaER", "CaOut", "RyRO"]:
        output_name = "all"
    elif  chosen_specie in ["STIM_CaER", "Orai"]:
        output_name = "RyR_Orai"

    try:
        specie_list = specie_dict[chosen_specie]
    except AttributeError:
        sys.exit("Unnkown specie %s" % s)
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04",
                "dend05", "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
    what_species = ["SOCE_", ""]
    dur_dict = {
        "SOCE_": " SOCE",
        "": " no SOCE"
    }
    
    organization = ["tubes"]
    fig1 = utils.make_distance_fig_det(directories, descr,
                                       dend_diam,
                                       stims, what_species,
                                       
                                       dur_dict,
                                       reg_list, output_name, 
                                       colors, types, symbol, marker)
    fig1.savefig("ER_distance_no_cam_det_lim.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("ER_distance_no_cam_det_lim.eps", dpi=100,
                 bbox_inches="tight")

