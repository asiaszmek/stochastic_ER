import os
import utility_functions as utils
from matplotlib.lines import Line2D


colors ={

        "1.2": 'tab:blue',
        "2.4": 'tab:purple',
        "6.0": 'tab:green'

}

directories = [
    ["Ca_wave_RyR2CaM_simple_SERCA_SOCE", "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_RyR2CaM_aging", "model_RyR2CaM_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_PMCA_aging", "model_aging%s_simple_SERCA_PMCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_normal_SERCA_aging", "model_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
]
legend_elements = [Line2D([0], [0], color='k', marker="o", fillstyle="full",
                          lw=0, label='ctrl'),
                   Line2D([0], [0], color="k", marker='s', fillstyle="none",
                          lw=0, label='old age with RyR2CaM'),
                   Line2D([0], [0], color="k", marker='s', fillstyle="full",
                          lw=0, label='ctrl + RyR2 no CaM'),
                   Line2D([0], [0], color='k', marker="o", fillstyle="none",
                          lw=0, label='old age'),]


stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]

if __name__ == '__main__':
    types = [ "ctrl",  "old age - RyR2 + RyR2CaM", "ctrl + RyR2",
              "old age",]

    marker = ["o", "s", "s", "o"]
    fillstyle = ["full", "none", "full", "none"]

    output_name = "all"
    fig1 = utils.make_distance_fig_sep_dends(directories,
                                             dend_diam,
                                             stims,
                                             output_name, 
                                             colors,
                                             types, marker, fillstyle,
                                             method="regular",
                                             legend=legend_elements)
    fig1.savefig("Aging_spatial_specificity_dissection.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("Aging_spatial_specificity_dissection.eps", dpi=100,
                 bbox_inches="tight")


