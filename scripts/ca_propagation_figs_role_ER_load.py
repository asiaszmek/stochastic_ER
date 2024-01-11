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

colors = ['tab:blue', 'tab:olive', 'tab:green', 'tab:red',  'tab:purple', 'tab:brown',
          'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', "b", "olive", "g",
          "r", "p" ]
stim_dend = "dend26"
directories = [
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER",
    "Ca_wave_simple_SERCA_no_SOCE_Breit_2018",
    "Ca_wave_simple_SERCA_no_SOCE_largerER",
]
descr = {
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE": "RyR2CaM",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER": "RyR2CaM_largerER",
    "Ca_wave_simple_SERCA_no_SOCE_Breit_2018": "RyR",
    "Ca_wave_simple_SERCA_no_SOCE_largerER": "RyR_largerER"
}
types = {
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE": "CaM",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER": "CaM + 110%ER",
    "Ca_wave_simple_SERCA_no_SOCE_Breit_2018": "no CaM",
    "Ca_wave_simple_SERCA_no_SOCE_largerER": "no CaM + 110%ER",
}
marker = {
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE": "full",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER": "full",
    "Ca_wave_simple_SERCA_no_SOCE_Breit_2018": "none",
    "Ca_wave_simple_SERCA_no_SOCE_largerER": "none",
}


symbol = {
    "tubes": "d",
    "baloon": "o",
}


dend_f = {
    "350 nM":
    [
        "model_%s_simple_SERCA_tubes_diam_1.2_um_50_um_0350_nM.h5",
        "model_%s_simple_SERCA_tubes_diam_2.4_um_50_um_0350_nM.h5",
        "model_%s_simple_SERCA_tubes_diam_6.0_um_50_um_0350_nM.h5",



           
    ],
       "700 nM":
    [
      
        "model_%s_simple_SERCA_tubes_diam_1.2_um_50_um_0700_nM.h5",        
        "model_%s_simple_SERCA_tubes_diam_2.4_um_50_um_0700_nM.h5",        
        "model_%s_simple_SERCA_tubes_diam_6.0_um_50_um_0700_nM.h5",        
    
       
    ],
       "1050 nM":
    [
        "model_%s_simple_SERCA_tubes_diam_1.2_um_50_um_1050_nM.h5",             
        "model_%s_simple_SERCA_tubes_diam_2.4_um_50_um_1050_nM.h5",        
        "model_%s_simple_SERCA_tubes_diam_6.0_um_50_um_1050_nM.h5",        
      
         
    ],
}
labels = {
    "350 nM":
    [
        "1.2 um",
        "2.4 um",        
        "6.0 um",

    ],
       "700 nM":
    [

        "1.2 um",
        "2.4 um",        
        "6.0 um",


    ],
       "1050 nM":
    [
        "1.2 um",
        "2.4 um",        
        "6.0 um",        

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
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
        
                          
    fig1 = utils.ca_wave_propagation_figs(directories, descr,
                                    dend_f, specie_list, reg_list,
                                    output_name, colors, labels,
                                    types, marker)

    fig1.savefig("ER_load_impact.eps",
                 bbox_inches="tight", pad_inches=0.2)
    fig1.savefig("ER_load_impact.png",
                 bbox_inches="tight", pad_inches=0.2)
  
