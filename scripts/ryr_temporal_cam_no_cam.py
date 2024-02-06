import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
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

directories = {
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE":"model_RyR2CaM%s_simple_SERCA%s_tubes_diam_%s_um_50_um_%s_nM.h5",
    "Ca_wave_simple_SERCA_SOCE":"model_RyR%s_simple_SERCA%s_tubes_diam_%s_um_50_um_%s_nM.h5"}




stims = [ "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]

    

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


if __name__ == '__main__':
    what_species = ["_SOCE", ""]
    organization = ["tubes"]
    types = {
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE": "ctrl",
    "Ca_wave_simple_SERCA_SOCE": "no CaM",
    }
    output_name = "all"
    fig1 = utils.make_decay_constant_fig(directories,
                                         dend_diam,
                                         stims, what_species,
                                         organization,
                                         dur_dict,
                                         output_name, 
                                         colors,
                                         types)
    fig1.savefig("CaM_no_CaM_temporal_short.png", dpi=100, bbox_inches="tight")
    fig1.savefig("CaM_no_CaM_temporal_short.eps", dpi=100, bbox_inches="tight")
    plt.show()
