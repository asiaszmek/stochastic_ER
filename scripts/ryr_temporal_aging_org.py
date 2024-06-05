import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors ={
    "Ca_wave_aging": {
        "1.2": 'tab:blue',
        "2.4": 'tab:purple',
        "6.0": 'tab:green'
    }
}

stim_dend = "dend26"


fname = "model_aging%s_simple_SERCA_%s_diam_%s_um_50_um_%s_nM.h5"


stims = [ "0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]
directory =  "Ca_wave_aging"


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

    organization = [ "tubes", "baloon"]

    types = {
    "baloon": "RyRCaM in EPJ",
    "tubes": "uniform RyRCaM",
    }
    
    output_name = "all"
    fig1 = utils.make_decay_constant_fig_ctrl(fname, directory,
                                              dend_diam,
                                              stims,
                                              organization,
                                              output_name, 
                                              colors,
                                              types)
    fig1.savefig("ER_org_temporal_aging_ER_org.png", dpi=100, bbox_inches="tight")
    fig1.savefig("ER_org_temporal_aging_ER_org.eps", dpi=100, bbox_inches="tight")
    plt.show()
