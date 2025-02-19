import os
import utility_functions as utils


colors ={

        "1.2": 'tab:blue',
        "2.4": 'tab:purple',
        "6.0": 'tab:green'

}

directories = [
    ["Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"], 
    ["Ca_wave_RyR2CaM_0.8_PMCA",
    "model_RyR2CaM%s_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_50_um_%s_nM.h5"] ,
    ["Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER",
    "model_RyR2CaM%s_largerER_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_RyR2CaM_120_PMCA",
    "model_RyR2CaM%s_simple_SERCA_SOCE_120_PMCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    
]
 
stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]

if __name__ == '__main__':
    types = [ "ctrl",  "80% PMCA kcat", "110% Ca in ER",
              "120% PMCA",]

    marker = ["o", "s", "s", "o"]
    fillstyle = ["full", "none", "full", "none"]

    output_name = "all"
    fig1 = utils.make_distance_fig_sep_dends(directories,
                                             dend_diam,
                                             stims,
                                             output_name, 
                                             colors,
                                             types, marker, fillstyle, method="regular")
    fig1.savefig("Param_sens_spatial_specificity_dissection.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("Param_sens_specificity_dissection.eps", dpi=100,
                 bbox_inches="tight")


