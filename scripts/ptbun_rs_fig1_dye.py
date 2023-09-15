import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils
colors_fura = ["k", "b", "g"]
colors_no_fura = ["tab:gray", "tab:blue", "tab:green"]


def set_ylim(ax, mini, maxi):
    for x in ax:
        x.set_ylim([mini + 0.01, maxi + 0.01])



def get_conc(fullname, specie_list, region_list, output_name):
   
    my_file = h5py.File(fullname)
    conc_dict = {}
    time_dict = {}
    
    for trial in my_file.keys():
        if trial == "model":
            continue
        conc, voxels = utils.get_dynamics_in_region(my_file,
                                                    specie_list,
                                                    region_list, trial,
                                                    output_name)
        conc_dict[trial] = conc
        time = utils.get_times(my_file, trial, output_name)
        time_dict[trial] = time
    lmin = min([len(conc) for conc in conc_dict.values()])
    time_end = min([time[-1] for time in time_dict.values()])
    time_len = min([len(time) for time in time_dict.values()])
    time = np.linspace(0, time_end, time_len)
    shape2 = max([conc.shape[1] for conc in conc_dict.values()])
    conc_mean = np.zeros((lmin, shape2))
    for conc in conc_dict.values():
        conc_mean[:lmin, :] += conc[:lmin, :]
    conc_mean /= len(conc_dict)
    return voxels, time, conc_mean
        


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
stims = ["0350", "0700", "1050"]
stim_labels = ["2 uM Ca injection", "4 uM Ca injection", "10 uM Ca injection"]
branch_diams = [1.2, 2.4, 6.0]
Fura_specie = "Fura2Ca"
t_stim = 3000  # sec
dye_base = "model_RyR_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nMFura2.h5"
nodye_base = "model_RyR_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(2, 3, figsize=(15, 12))

#  Dye figs
mini = 200000
maxi = 0 
for i, x in enumerate(ax[0]):
    x.set_title(stim_labels[i], fontsize=20)
    for k, b_diam in enumerate(branch_diams):
        fname = dye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, basic_RyR_no_SOCE_dir, fname)
        voxels, time, conc_mean = get_conc(full_name, ["Fura2Ca"], reg_list,
                                           output_name)
        mean = conc_mean[:, :int(t_stim/(time[1]-time[0]))].mean(axis=0)
        conc_mean = (conc_mean - mean)/mean
        # mean fluo
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(conc_mean.shape[1]):
            max_fluo_vals[j] = conc_mean[:, j].max()
        mini = min(mini, min(max_fluo_vals))
        maxi = max(maxi, max(max_fluo_vals))
        x.plot(vox_axis, max_fluo_vals, colors_fura[k],
               label="dend diam=%2.1f um" % b_diam)
    if not i:
        x.legend()
        x.set_ylabel("Fluorescence", fontsize=20)
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
        full_name = os.path.join(cur_dir, basic_RyR_no_SOCE_dir, fname)
        voxels, time, ca = get_conc(full_name, ["Ca"], reg_list, output_name)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(conc_mean.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x.plot(vox_axis, max_fluo_vals, colors_fura[k],
               label="%2.1f um + Fura2" % b_diam)
        
        maxi = max(maxi, max(max_fluo_vals))
        mini = min(mini, min(max_fluo_vals))

        fname_no_dye = nodye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, basic_RyR_no_SOCE_dir, fname_no_dye)
        voxels, time, ca = get_conc(full_name, ["Ca"], reg_list, output_name)
        vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
        max_fluo_vals = np.zeros_like(vox_axis)
        for j in range(conc_mean.shape[1]):
            max_fluo_vals[j] = ca[:, j].max()/1000
        x.plot(vox_axis, max_fluo_vals, colors_no_fura[k],
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
fig.savefig("Ca_dye_effects.png", dpi=300, bbox_inches="tight")
plt.show()
