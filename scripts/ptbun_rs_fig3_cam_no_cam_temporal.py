import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {"no_SOCE_no_CaM": "tab:blue",
          "_no_CaM": "tab:green",
          "no_SOCE_CaM": "tab:cyan",
          "_CaM": "tab:olive",
}
labels = {"no_SOCE_no_CaM": "no SOCE RyR dis-inh.",
          "_no_CaM": "RyR dis-inh.",
          "no_SOCE_CaM": "no SOCE",
          "_CaM": "ctrl",
}

def set_ylim(ax, mini, maxi):
    for x in ax:
        x.set_ylim([mini + 0.01, maxi + 0.05])


reg_list = ["dend25","dend26", "dend27"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
basic_RyR_no_SOCE_dir = "Ca_wave_simple_SERCA_no_SOCE_Breit_2018"
basic_RyR_SOCE_dir = "Ca_wave_simple_SERCA_SOCE"
RyRCaM_no_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE"
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stims = ["0350", "0700"]
stim_labels = ["2 uM Ca injection", "4 uM Ca injection", "10 uM Ca injection"]
branch_diams = [1.2, 2.4, 6.0]
Fura_specie = "Fura2Ca"
t_start = 3000
idx_start = t_start
base_noSOCE = "model_RyR_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
base_SOCE = "model_RyR_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

CaM_noSOCE = "model_RyR2CaM_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))

mini = []
maxi = []


for i, x in enumerate(ax):
   
    for k, b_diam in enumerate(branch_diams):
        x[k].set_title(stim_labels[i]+ " diam %2.1f um" % b_diam,
                       fontsize=12)
        fname = base_noSOCE % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, basic_RyR_no_SOCE_dir, fname)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        output = ca.mean(axis=1)
        x[k].plot(time, output, colors["no_SOCE_no_CaM"],
                  label=labels["no_SOCE_no_CaM"])
        maxi.append(max(output))
        mini.append(min(output))

        fname_SOCE = base_SOCE % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, basic_RyR_SOCE_dir, fname_SOCE)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        output = ca.mean(axis=1)
        x[k].plot(time, output, colors["_no_CaM"],
                  label=labels["_no_CaM"])
        maxi.append(max(output))
        mini.append(min(output))

        
        x[k].tick_params(labelsize=14)    
        fname = CaM_noSOCE % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_no_SOCE_dir, fname)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        output = ca.mean(axis=1)
        x[k].plot(time, output, colors["no_SOCE_CaM"],
                  label=labels["no_SOCE_CaM"])
        maxi.append(max(output))
        mini.append(min(output))

        
        fname_SOCE = CaM_SOCE % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname_SOCE)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        output = ca.mean(axis=1)
        x[k].plot(time, output, colors["_CaM"],
                  label=labels["_CaM"])
        maxi.append(max(output))
        mini.append(min(output))

        if not k:
            x[k].legend()
            x[k].set_ylabel("Calcium [uM]", fontsize=20)
        else:
            x[k].set_yticks([])

        if i < 2:
            x[k].set_xticks([])
        else:
            x[k].set_xlabel("Distance from stim [um]", fontsize=20)
            x[k].tick_params(labelsize=14)


#ax[2, 2].legend()
set_ylim(ax[1], min(mini), max(maxi))
set_ylim(ax[0], min(mini), max(maxi))
#set_ylim(ax[2], min(mini), max(maxi))

fig.savefig("Ca_decay_stim_point.svg", dpi=100, bbox_inches="tight")
plt.show()
