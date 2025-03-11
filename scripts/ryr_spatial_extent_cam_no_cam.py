import os
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import utility_functions as utils

import matplotlib.legend as mlegend
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

plt.rcParams['text.usetex'] = True

colors = {1.2: "tab:blue",
          2.4: "tab:purple",
          6.0: "tab:green",
}
labels = {"no_SOCE_no_CaM": "no CaM no SOCE",
          "_no_CaM": "no CaM",
          "no_SOCE_CaM": "ctrl no SOCE",
          "_CaM": "ctrl",
}


reg_list = ["dend25", "dend26", "dend27"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
basic_RyR_SOCE_dir = "Ca_wave_simple_SERCA_SOCE"
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stims = ["0175", "0350", "0700", "1050", "2000"]
stim_label = "2 uM Ca injection"
branch_diams = [1.2, 2.4, 6.0]

injections = {
    1.2:{
        "0175":{
            "":20*1200,
            "_3s_injection":1.5*12000,
        },

        "0350":{
            "":40*1200,
            "_3s_injection":3*12000,
        },
        "0700":{
            "":40*2400,
            "_3s_injection":3*24000,
        },
        "1050":{
            "":40*4000,
            "_3s_injection":3*40000,
        },
         "2000":{
            "":40*4800,
            "_3s_injection":3*48000,
        },
    },
    2.4:{
        "0175":{
            "":20*2000,
            "_3s_injection":1.5*20000,
        },

        "0350":{
            "":40*2000,
            "_3s_injection":3*20000,
        },
        "0700":{
            "":40*4000,
            "_3s_injection":3*40000,
        },
        "1050":{
            "":40*6000,
            "_3s_injection":3*60000,
        },
        "2000":{
            "":40*8000,
            "_3s_injection":3*80000,
        },
    },
    6.0:{
        "0175":{
            "":20*4000,
            "_3s_injection":1.5*40000,
        },

        "0350":{
            "":40*4000,
            "_3s_injection":3*40000,
        },
        
        "0700":{
            "": 40*8000,
            "_3s_injection":3*80000,
        },
        "1050":{
            "":40*12000,
            "_3s_injection":3*120000,
        },
        "2000":{
            "":40*16000,
            "_3s_injection":3*160000,
        },
    },
}

t_start = 3000
idx_start = t_start

base_SOCE = "model_RyR%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fname = {
    "_CaM": CaM_SOCE,
    "_no_CaM": base_SOCE
    }

directory ={
    "_CaM": "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "_no_CaM": "Ca_wave_simple_SERCA_SOCE",
    }

stim_types = [ ""]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

suffix = {"": "",
          #"_3s_injection": " bAP"
          }


for k, b_diam in enumerate(branch_diams):
    for inh in ["_no_CaM", "_CaM"]:
        for j, s in enumerate(stim_types):
            y = []
            y_err = []
            x = []
            for stim in stims:
                fname_SOCE = fname[inh] % (s, b_diam, stim)
                full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
                ca_dict, time_dict = utils.get_conc(full_name, ["Ca"],
                                                    reg_list,
                                                    output_name)
                try:
                    ca_out = utils.get_array(ca_dict, "Ca")
                except ValueError:
                    continue
                ca = np.array(ca_out)/1000
                ca = ca.mean(axis=1)
                output = ca.max(axis=1)
                y.append(np.mean(output))
                x.append(injections[b_diam][stim][s]/b_diam)
                y_err.append(output.std()/len(output)**0.5)
            print(x, y, y_err)
            if inh == "_CaM" and s == "":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="d",
                            linestyle="")
            elif inh == "_no_CaM" and s == "":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="d",
                            linestyle="",
                            fillstyle="none")

ims1 = Line2D([0], [0], color="tab:blue", marker="d", fillstyle="full",
              lw=0)
ims2 = Line2D([0], [0], color="tab:purple", marker="d", fillstyle="full",
              lw=0)
ims3 = Line2D([0], [0], color="tab:green", marker="d", fillstyle="full",
              lw=0)

ims4 = Line2D([0], [0], color="tab:blue", marker="d", fillstyle="none",
              lw=0)
ims5 = Line2D([0], [0], color="tab:purple", marker="d", fillstyle="none",
              lw=0)
ims6 = Line2D([0], [0], color="tab:green", marker="d", fillstyle="none",
              lw=0)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.set_ylabel(r"max(Ca) $(\mathrm{\mu M})$", fontsize=15)
ax.set_xlabel(r"total injected ions/diam $(\mathrm{1/\mu m})$", fontsize=15)
extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
legend_handle = [extra, extra, extra, extra, ims1, ims4, extra, ims2, ims5, extra, ims3, ims6]
label_col_1 = ["diam"]
label_j_1 = ["ctrl"]
label_j_2 = ["no CaM"]

label_empty = [""]
legend_labels = np.concatenate([label_col_1, label_j_1, label_j_2,
                                [r"1.2 $\mathrm{\mu m}$"], label_empty * 2,
                                [r"2.4 $\mathrm{\mu m}$"], label_empty * 2,
                                [r"6.0 $\mathrm{\mu m}$"], label_empty * 2])
ax.legend(legend_handle, legend_labels, 
          ncol = 4, shadow = True, handletextpad = -2)

#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.tick_params(axis='both', which='major', labelsize=14)

fig.savefig("Spatial_sum_cam.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Spatial_sum_cam.png", dpi=100,
            bbox_inches="tight")


