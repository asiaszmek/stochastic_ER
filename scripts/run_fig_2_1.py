import os
import utility_functions as utils
from matplotlib.lines import Line2D

legend_elements = [Line2D([0], [0], color='k', marker="o", fillstyle="full",
                          lw=0, label='ctrl'),
                   Line2D([0], [0], color="k", marker='o', fillstyle="none",
                          lw=0, label="no SOCE"),
                  ]

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


types = [ "RyR2CaM+SOCE", "RyR2CaM"]
marker = ["o", "o"]
fillstyle = ["full", "none"]
stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["2.4"] #["1.2", "2.4", "6.0"]


if __name__ == '__main__':
    fig1 = utils.make_distance_fig_sep_dends(directories, dend_diam,
                                             stims, "all", 
                                             colors, types, marker,
                                             fillstyle, legend=legend_elements)
    fig1.savefig("ER_distance_SOCE_no_SOCE.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("ER_distance_SOCE_no_SOCE.eps", dpi=100,
                 bbox_inches="tight")
