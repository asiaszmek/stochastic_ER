import os
import utility_functions as utils
from matplotlib.lines import Line2D

colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}
directories = [
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
        "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE_baloon",
        "model_RyR2CaM%s_simple_SERCA_SOCE_baloon_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_normal_SERCA_aging",
        "model_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_normal_SERCA_aging_baloon",
        "model_aging%s_simple_SERCA_baloon_diam_%s_um_50_um_%s_nM.h5"
    ]
]

stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]

legend_elements = [Line2D([0], [0], color='k', marker="o", fillstyle="full",
                          lw=0, label='uniform RyR2CaM'),
                   Line2D([0], [0], color="k", marker='^', fillstyle="full",
                          lw=0, label="RyR2CaM in EPJ"),
                   Line2D([0], [0], color="k", marker='o', fillstyle="none",
                          lw=0, label="old age uniform"),
                   Line2D([0], [0], color='k', marker="^", fillstyle="none",
                          lw=0, label='old age EPJ'),]



if __name__ == '__main__':
    types = ["uniform RyR2CaM", "RyR2CaM in EPJ", "old age uniform RyR2 and RyR2CaM", "old age RyR2 and RyR2CaM in EPJ"]
    markers = ["o", "^", "o", "^"]
    fillstyle = ["full", "full", "none", "none"]

    fig1 = utils.make_distance_fig_sep_dends(directories,
                                             dend_diam,
                                             stims,
                                             "all", 
                                             colors,
                                             types, marker=markers,
                                             fillstyle=fillstyle,
                                             legend=legend_elements)
    fig1.savefig("Aging_distance.png", dpi=100, bbox_inches="tight")
    fig1.savefig("Aging_distance.eps", dpi=100, bbox_inches="tight")
