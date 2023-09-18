import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {"_no_CaM": "tab:green",
          "_CaM": "tab:olive",
          "_no_CaM_ove": "tab:pink",
          "_CaM_ove": "tab:purple",
}
labels = {"_no_CaM": "RyR dis-inh.",
          "_CaM": "ctrl",
          "_no_CaM_ove": "RyR over-exp. and dis-inh.",
          "_CaM_ove": "RyR over-exp.",
}

def set_ylim(ax, mini, maxi):
    for x in ax:
        x.set_ylim([mini + 0.01, maxi + 0.05])
   


reg_list = ["dend", "dend01", "dend02", "dend03", "dend04", "dend05",
            "dend06", "dend07", "dend08", "dend09", "dend10",
            "dend11", "dend12", "dend13", "dend14", "dend15",
            "dend16", "dend17", "dend18", "dend19", "dend20",
            "dend21", "dend22", "dend23", "dend24", "dend25",
            "dend26", "dend27", "dend28", "dend29", "dend30",
            "dend31", "dend32", "dend33", "dend34", "dend35",
            "dend36", "dend37", "dend38", "dend39", "dend40",
            "dend41", "dend42", "dend43", "dend44", "dend45",
            "dend46", "dend47", "dend48", "dend49", "dend50",
            "dend51"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
basic_RyR_no_SOCE_dir = "Ca_wave_simple_SERCA_no_SOCE_Breit_2018"
basic_RyR_SOCE_dir = "Ca_wave_simple_SERCA_SOCE"
RyRCaM_no_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE"
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stims = ["0350"]
stim_labels = ["2 uM Ca injection", "4 uM Ca injection", "10 uM Ca injection"]
branch_diams = [2.4, 6.0]
Fura_specie = "Fura2Ca"
ryr_exp = ["", "_nc"]
ryr_exp_label = ["", "RyR overexpression"]
base_noSOCE = "model_RyR_simple_SERCA%s_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
base_SOCE = "model_RyR_simple_SERCA_SOCE%s_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

CaM_noSOCE = "model_RyR2CaM_simple_SERCA%s_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM_simple_SERCA_SOCE%s_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))

mini = []
maxi = []


for i, x in enumerate([ax]):
   
    for k, b_diam in enumerate(branch_diams):
    
        x[k].set_title(stim_labels[i]+ " diam %2.1f um" % b_diam,
                       fontsize=12)
        fname = base_SOCE % (ryr_exp[0], b_diam, stims[i])
        full_name = os.path.join(cur_dir, basic_RyR_SOCE_dir, fname)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(ca.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x[k].plot(vox_axis, max_fluo_vals, colors["_no_CaM"],
                  label=labels["_no_CaM"])
        maxi.append(max(max_fluo_vals))
        mini.append(min(max_fluo_vals))

        fname_SOCE = base_SOCE % (ryr_exp[1], b_diam, stims[i])
        full_name = os.path.join(cur_dir, basic_RyR_SOCE_dir, fname_SOCE)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(ca.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x[k].plot(vox_axis, max_fluo_vals, colors["_no_CaM_ove"],
                  label=labels["_no_CaM_ove"])
        maxi.append(max(max_fluo_vals))
        mini.append(min(max_fluo_vals))

        
        x[k].tick_params(labelsize=14)    
        fname = CaM_SOCE % (ryr_exp[0], b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(ca.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x[k].plot(vox_axis, max_fluo_vals, colors["_CaM"],
                  label=labels["_CaM"])
        maxi.append(max(max_fluo_vals))
        mini.append(min(max_fluo_vals))

        
        fname_SOCE = CaM_SOCE % (ryr_exp[1], b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname_SOCE)
        voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(ca.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x[k].plot(vox_axis, max_fluo_vals, colors["_CaM_ove"],
                  label=labels["_CaM_ove"])
        maxi.append(max(max_fluo_vals))
        mini.append(min(max_fluo_vals))

        if not k:
            x[k].legend()
            x[k].set_ylabel("Calcium [uM]", fontsize=20)
        else:
            x[k].set_yticks([])

        # if i < 2:
        #     x[k].set_xticks([])
        # else:
        x[k].set_xlabel("Distance from stim [um]", fontsize=20)
        x[k].tick_params(labelsize=14)


ax[1].legend()
set_ylim(ax, min(mini), max(maxi))
# set_ylim(ax[0], min(mini), max(maxi))
# set_ylim(ax[2], min(mini), max(maxi))

fig.savefig("overexpression_effects.svg", dpi=100, bbox_inches="tight")
plt.show()
