import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils

colors = ["tab:blue", "tab:orange", "tab:green"]


def set_ylim(ax, mini, maxi):
    for x in ax:
        x.set_ylim([mini + 0.01, maxi + 0.01])


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
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stims = ["0350", "0700", "1050"]
stim_labels = ["2 uM Ca injection", "4 uM Ca injection", "10 uM Ca injection"]
branch_diams = [1.2, 2.4, 6.0]
Fura_specie = "Fura2Ca"
t_stim = 3000  # sec
dye_base = "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.xmlFura2.h5"
nodye_base = "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(2, 3, figsize=(15, 12))

#  Dye figs
mini = 200000
maxi = 0 
for i, x in enumerate(ax[0]):
    x.set_title(stim_labels[i], fontsize=20)
    for k, b_diam in enumerate(branch_diams):
        fname = dye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname)
        try:
            voxels, time, conc_mean = utils.get_conc(full_name, ["Fura2Ca"],
                                                     reg_list, output_name)
        except FileNotFoundError or OSError:
            continue
        mean = conc_mean[:, :int(t_stim/(time[1]-time[0]))].mean(axis=0)
        conc_mean = (conc_mean - mean)/mean
        # mean fluo
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(conc_mean.shape[1]):
            max_fluo_vals[j] = conc_mean[:, j].max()
        mini = min(mini, min(max_fluo_vals))
        maxi = max(maxi, max(max_fluo_vals))
        x.plot(vox_axis, max_fluo_vals, colors[k], marker="d",
               linestyle="", fillstyle="none",
               label="dend diam=%2.1f um" % b_diam)
    if not i:
        x.legend()
        x.set_ylabel("%Fluorescence", fontsize=20)
    else:
        x.set_yticklabels([])
    x.set_xticklabels([])
    x.tick_params(labelsize=14)

set_ylim(ax[0], mini, maxi)
mini = 200000
maxi = 0 
for i, x in enumerate(ax[1]):

    for k, b_diam in enumerate(branch_diams):
        fname = dye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname)
        try:
            voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        except FileNotFoundError or OSError:
            continue
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(conc_mean.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x.plot(vox_axis, max_fluo_vals, colors[k], marker="d",
               linestyle="", fillstyle="none",
               label="%2.1f um + Fura2" % b_diam)
        
        maxi = max(maxi, max(max_fluo_vals))
        mini = min(mini, min(max_fluo_vals))
        
        fname_no_dye = nodye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname_no_dye)
        try:
            voxels, time, ca = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
        except FileNotFoundError or OSError:
            continue
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(conc_mean.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x.plot(vox_axis, max_fluo_vals, colors[k],
                marker="d",
               linestyle="", fillstyle="full",
               label="%2.1f um - Fura2" % b_diam)
        maxi = max(maxi, max(max_fluo_vals))
        mini = min(mini, min(max_fluo_vals))

        
    
    if not i:
        x.legend()
        x.set_ylabel("Calcium [uM]", fontsize=20)
    else:
        x.set_yticklabels([])
    x.tick_params(labelsize=14)
    x.set_xlabel("Distance from stim [um]", fontsize=20)







set_ylim(ax[1], mini, maxi)
fig.savefig("Ca_dye_effects.eps", dpi=100, bbox_inches="tight")
fig.savefig("Ca_dye_effects.png", dpi=100, bbox_inches="tight")
plt.show()
