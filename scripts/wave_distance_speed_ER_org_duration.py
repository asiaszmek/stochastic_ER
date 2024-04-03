import os
import h5py
import numpy as np
import sys
import argparse
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
import utility_functions as utils

marker = {
    "": "d",
    "_3s_injection": "o"
}

dur_dict = {
    "": " 40 ms",
    "_3s_injection": " 3 ms"
}
              

colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}

stim_dend = "dend26"

directories = [
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
]

descr = {
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE": "RyR2CaM",
}

types = {
    "baloon": "RyRCaM in dendritic membrane",
    "tubes": "uniformly distributed RyRCaM",
}

stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]

    
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


    fname = "model_%s%s_simple_SERCA_SOCE_%s_diam_%s_um_50_um_%s_nM.h5"
    what_species = ["", "_3s_injection"]
    organization = ["baloon", "tubes"]
    fig1 = utils.make_distance_fig_2_4(fname, directories, descr,
                                       dend_diam,
                                       stims, what_species,organization,
                                       dur_dict,
                                       reg_list, output_name, 
                                       colors, types, marker)
    fig1.savefig("ER_org_distance.png", dpi=100, bbox_inches="tight")
    fig1.savefig("ER_org_distance.eps", dpi=100, bbox_inches="tight")
