import os
import sys
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import utility_functions as utils
from scipy.constants import Avogadro


t_init = 3000

colors = ['tab:blue', 'tab:olive', 'tab:green', 'tab:blue', 'tab:olive', 'tab:green']
stim_dend = "dend26"
directories = [
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    #"Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
]

symbol = {
    "tubes": "d",
    "baloon": "o",
}
types = {"Ca_wave_RyR2CaM_simple_SERCA_SOCE": ""}
marker = {"Ca_wave_RyR2CaM_simple_SERCA_SOCE": "full"}
dend_f = {
    "350 nM":
    [
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_1.2_um_50_um_0350_nM.h5",
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_0350_nM.h5",
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_6.0_um_50_um_0350_nM.h5",
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_1.2_um_50_um_0350_nM.h5",
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_2.4_um_50_um_0350_nM.h5",
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_6.0_um_50_um_0350_nM.h5",


           
    ],
       "700 nM":
    [
      
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_1.2_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_6.0_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_1.2_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_2.4_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_6.0_um_50_um_0700_nM.h5",        
    
       
    ],
       "1050 nM":
    [
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_1.2_um_50_um_1050_nM.h5",               
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_1050_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_tubes_diam_6.0_um_50_um_1050_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_1.2_um_50_um_1050_nM.h5",               
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_2.4_um_50_um_1050_nM.h5",        
        "model_RyR2CaM_3s_injection_simple_SERCA_SOCE_baloon_diam_6.0_um_50_um_1050_nM.h5",        
      
         
    ],
}
labels = {
    "350 nM":
    [
        "RyR2 uniform  1.2 um",
        "RyR2 uniform  2.4 um",        
        "RyR2 uniform  6.0 um",
        "RyR2 dendritic membrane  1.2 um",
        "RyR2 dendritic membrane  2.4 um",        
        "RyR2 dendritic membrane  6.0 um",        

    ],
       "700 nM":
    [

        "RyR2 uniform  1.2 um",
        "RyR2 uniform  2.4 um",        
        "RyR2 uniform  6.0 um",
        "RyR2 dendritic membrane  1.2 um",
        "RyR2 dendritic membrane  2.4 um",        
        "RyR2 dendritic membrane  6.0 um",        


    ],
       "1050 nM":
    [
        "RyR2 uniform  1.2 um",
        "RyR2 uniform  2.4 um",        
        "RyR2 uniform  6.0 um",        
        "RyR2 dendritic membrane  1.2 um",
        "RyR2 dendritic membrane  2.4 um",        
        "RyR2 dendritic membrane  6.0 um",        

    ],
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
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
    
    fig1 = utils.ca_wave_propagation_figs_bal_tubes(directories, "",
                                              dend_f, specie_list,
                                              reg_list, output_name,
                                              colors, labels, types, marker)
    fig1.savefig("Ca_wave_vs_distance_3_ms_injection.png",
                 bbox_inches="tight", pad_inches=0.1)
    fig1.savefig("Ca_wave_vs_distance_3_ms_injection.eps",
                 bbox_inches="tight", pad_inches=0.1)       
