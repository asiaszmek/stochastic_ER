import os
import utility_functions as utils

colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}


directories = [
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
        "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
        "model_RyR2CaM%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
]


types = [ "ctrl", "ctrl no SOCE"]
marker = ["o", "o"]
fillstyle = ["full", "none"]
stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]


if __name__ == '__main__':
    fig1 = utils.make_injection_vs_min_CaER(directories, dend_diam,
                                            stims, "all", 
                                            colors, types, marker,
                                            fillstyle)
    fig1.savefig("min_CaER_SOCE_no_SOCE.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("min_CaER_SOCE_no_SOCE.eps", dpi=100,
                 bbox_inches="tight")
