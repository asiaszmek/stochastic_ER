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
    ["Ca_wave_simple_SERCA_no_SOCE_Breit_2018", "model_RyR%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_RyR2CaM_simple_SERCA_SOCE_low_PMCA", "model_RyR2CaM%s_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_RyR2CaM_simple_SERCA_SOCE_2x_buffers", "model_RyR2CaM%s_simple_SERCA_SOCE_2x_buffer_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_RyR2CaM_aging", "model_RyR2CaM_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_normal_SERCA_aging", "model_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
]
legend_elements = [Line2D([0], [0], color='k', marker="o", fillstyle="full",
                          lw=0, label='ctrl'),
                   Line2D([0], [0], color="k", marker='o', fillstyle="none",
                          lw=0, label='ctrl - RyR2CaM + RyR2 no CaM'),
                   Line2D([0], [0], color="k", marker='^', fillstyle="full",
                          lw=0, label='ctrl with 80\% PMCA activity'),
                   Line2D([0], [0], color="k", marker='^', fillstyle="none",
                          lw=0, label='ctrl with 200\% Ca-buffers'),
                   Line2D([0], [0], color="k", marker='s', fillstyle="full",
                          lw=0, label='old age - RyR2 no CaM + RyR2CaM'),
                   Line2D([0], [0], color='k', marker="s", fillstyle="none",
                          lw=0, label='old age'),]


stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2"]#, "2.4", "6.0"]

if __name__ == '__main__':
    types = [ "ctrl",  "ctrl - RyR2CaM + RyR2 no CaM", "ctrl with low PMCA activity" ,"ctrl with 200\% Ca buffers","old age - RyR2 no CaM + RyR2CaM",
              "old age",]
    marker = ["o", "o", "^", "^", "s","s", ]
    fillstyle = ["full", "none", "full", "none", "full", "none"]

    output_name = "all"
    fig1 = utils.make_decay_constant_fig_sep_dends(directories,
                                                   dend_diam,
                                                   stims,
                                                   output_name, 
                                                   colors,
                                                   types, marker,
                                                   fillstyle)#,
                                                   #legend=legend_elements)
    fig1.savefig("Aging_temporal_dissection.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("Aging_temporal_dissection.eps", dpi=100,
                 bbox_inches="tight")


