import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils

colors_ctrl = ["tab:blue", "tab:purple", "tab:green"]
colors_dye = ["tab:cyan", "tab:pink", "tab:olive"]

def get_ylim(ax):
    mini = min([min(x.get_ylim()) for x in ax])
    maxi = max([max(x.get_ylim()) for x in ax])
    return mini, maxi


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
dye_base = "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_%2.1f_um_50_um_%s_nM.xmlFura2.h5"
nodye_base = "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(2, 3, figsize=(15, 12))

for i, x in enumerate(ax[0]):
    x.set_title(stim_labels[i], fontsize=20)
    for k, b_diam in enumerate(branch_diams):
        fname = dye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname)
        my_file = h5py.File(full_name)
        my_grid = utils.get_grid_list(my_file)
        vox_ind, vols = utils.get_dend_indices(my_grid, region=reg_list)
        voxels = sorted(vox_ind.keys())
        conc_dict, time_dict = utils.get_conc(full_name, ["Fura2Ca", "Ca"],
                                              reg_list,
                                              output_name)
        dt = time_dict["trial0"][1] - time_dict["trial0"][0]
        ca_out = [conc_dict["Ca"][key]
                  for key in conc_dict["Ca"].keys()]
        ca = np.array(ca_out)/1000
        fura2 = [conc_dict["Fura2Ca"][key]
                 for key in conc_dict["Fura2Ca"].keys()]
        fura2 = np.array(fura2)
        
        mean = fura2[:, :int(3000/dt)].mean(axis=2)
        fura_scaled = (fura2 - mean[:, :, None])/mean[:, :, None]
        fluo_vals = fura_scaled.max(axis=2)
        max_fluo_vals = fluo_vals.mean(axis=0)
        max_fluo_error = fluo_vals.std(axis=0)/fluo_vals.shape[0]**0.5
        
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        x.plot(vox_axis, max_fluo_vals, color=colors_dye[k], 
               label="dend diam=%2.1f um" % b_diam)
        x.fill_between(vox_axis, max_fluo_vals-max_fluo_error,
                    max_fluo_vals+max_fluo_error,
                    color=colors_dye[k], alpha=0.2)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        ca_vals = ca.max(axis=2)
        max_ca_vals = ca_vals.mean(axis=0)
        max_ca_error = ca_vals.std(axis=0)/len(ca_vals)**0.5
        ax[1][i].plot(vox_axis, max_ca_vals, color=colors_dye[k],
                      label="%2.1f um + Fura2" % b_diam)
        ax[1][i].fill_between(vox_axis, max_ca_vals-max_ca_error,
                              max_ca_vals+max_ca_error, color=colors_ctrl[k],
                              alpha=0.2)

    if not i:
        x.legend()
        x.set_ylabel("%Fluorescence", fontsize=20)
    else:
        x.set_yticklabels([])
    x.set_xticklabels([])
    x.tick_params(labelsize=14)

mini, maxi = get_ylim(ax[0])
set_ylim(ax[0], mini, maxi)

for i, x in enumerate(ax[1]):
    for k, b_diam in enumerate(branch_diams):
        fname_no_dye = nodye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname_no_dye)
        my_file = h5py.File(full_name)
        my_grid = utils.get_grid_list(my_file)
        vox_ind, vols = utils.get_dend_indices(my_grid, region=reg_list)
        voxels = sorted(vox_ind.keys())
        conc_dict, time_dict = utils.get_conc(full_name, ["Ca"],
                                              reg_list,
                                              output_name)
        dt = time_dict["trial0"][1] - time_dict["trial0"][0]
        ca_out = [conc_dict["Ca"][key]
                  for key in conc_dict["Ca"].keys()]
        ca = np.array(ca_out)/1000

        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        ca_vals = ca.max(axis=2)
        max_ca_vals = ca_vals.mean(axis=0)
        max_ca_error = ca_vals.std(axis=0)/len(ca_vals)**0.5
        x.plot(vox_axis, max_ca_vals,
               color=colors_ctrl[k],
               label="%2.1f um - Fura2" % b_diam)
        x.fill_between(vox_axis, max_ca_vals-max_ca_error,
                       max_ca_vals+max_ca_error, color=colors_ctrl[k],
                       alpha=0.2)
    
    if not i:
        x.legend()
        x.set_ylabel("Calcium [uM]", fontsize=20)
    else:
        x.set_yticklabels([])
    x.tick_params(labelsize=14)
    x.set_xlabel("Distance from stim [um]", fontsize=20)
mini, maxi = get_ylim(ax[1])
set_ylim(ax[1], mini, maxi)



#  Dye figs
fig.savefig("Ca_dye_effects_baloon.eps", dpi=100, bbox_inches="tight")
fig.savefig("Ca_dye_effects_baloon.png", dpi=100, bbox_inches="tight")

