import os
import utility_functions as utils


colors = {"Ca_wave_RyR2CaM_simple_SERCA_SOCE": 'k',
          "Ca_wave_simple_SERCA_SOCE":"orange",
          "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE":"grey",
          "Ca_wave_simple_SERCA_no_SOCE_Breit_2018":"yellow",
}
stim_dend = "dend26"

directories = {
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE":"model_RyR2CaM%s_simple_SERCA%s_tubes_diam_%s_um_50_um_%s_nM.h5",
    
    "Ca_wave_simple_SERCA_SOCE":"model_RyR%s_simple_SERCA%s_tubes_diam_%s_um_50_um_%s_nM.h5",
    
}
stims = ["0175", "0350", "0700", "1050", "2000"]
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
        "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE": "ctrl no SOCE",
        "Ca_wave_simple_SERCA_no_SOCE_Breit_2018":"no CaM no SOCE",
    }
    output_name = "all"
    fig1 = utils.make_spatiotemporal_specificity_fig_sep_dends(directories,
                                                               dend_diam,
                                                               stims,
                                                               what_species,
                                                               output_name, 
                                                               colors,
                                                               types,
                                                               method="regular")
    fig1.savefig("CaM_spatiotemporal_dissection.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("CaM_spatiotemporal_dissection.eps", dpi=100,
                 bbox_inches="tight")
    fig1 = utils.make_spatiotemporal_specificity_fig_sep_dends(directories,
                                                               dend_diam,
                                                               stims,
                                                               what_species,
                                                               output_name, 
                                                               colors,
                                                               types,
                                                               method="fitexp")
    fig1.savefig("CaM_spatiotemporal_dissection_fitexp.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("CaM_spatiotemporal_dissection_fitexp.eps", dpi=100,
                 bbox_inches="tight")
