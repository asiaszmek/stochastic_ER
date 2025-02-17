import os
import utility_functions as utils


colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}
dir_list = [
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "Ca_wave_RyR2CaM_0.8_PMCA",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER",
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "Ca_wave_RyR2CaM_120_PMCA",
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "Ca_wave_2xRyR2CaM_simple_SERCA_SOCE",
]    
directories = {
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE":
    "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5", 
    "Ca_wave_RyR2CaM_0.8_PMCA":
    "model_RyR2CaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_50_um_%s_nM.h5" ,
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE":
    "model_RyR2CaM_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER":
    "model_RyR2CaM_largerER_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5",
    "Ca_wave_RyR2CaM_120_PMCA":
    "model_RyR2CaM_simple_SERCA_SOCE_120_PMCA_tubes_diam_%s_um_50_um_%s_nM.h5",
    "Ca_wave_2xRyR2CaM_simple_SERCA_SOCE":
    "model_2xRyR2CaM_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5", 
}

types = {
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE": "CaM",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER": "CaM + 110%ER",
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE": "CaM",
    "Ca_wave_RyR2CaM_0.8_PMCA": "CaM + 80% PMCA activity",
    "Ca_wave_RyR2CaM_120_PMCA": "CaM + 120% PMCA",
     "Ca_wave_2xRyR2CaM_simple_SERCA_SOCE": "CaM + 200% RyR2",
}


stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]


if __name__ == '__main__':
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend0", "dend04",
                "dend05", "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
    output_name = "all"        
    fig1 = utils.make_distance_fig_ratio(dir_list, directories, dend_diam,
                                         stims, ["Ca"], reg_list, output_name,
                                         colors, types)
    fig1.savefig("Param_var_distance_tubes.png", dpi=100, bbox_inches="tight")
    fig1.savefig("Param_var_distance_tubes.eps", dpi=100, bbox_inches="tight")

