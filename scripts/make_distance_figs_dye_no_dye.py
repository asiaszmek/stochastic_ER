import os
import h5py
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils

colors_dye = ["tab:cyan", "tab:pink", "tab:olive"]
colors_ctrl = ["tab:blue", "tab:purple", "tab:green"]

def get_fluo(full_name):
    try:
        my_file = h5py.File(full_name)
    except FileNotFoundError:
        return
    my_grid = utils.get_grid_list(my_file)
    vox_ind, vols = utils.get_dend_indices(my_grid, region=reg_list)
    voxels = sorted(vox_ind.keys())
    conc_dict, time_dict = utils.get_conc(full_name, ["Fura2Ca"],
                                          reg_list,
                                          output_name)
    try:
        dt = time_dict["trial0"][1] - time_dict["trial0"][0]
    except IndexError:
        return
    fura2 = utils.get_array(conc_dict,  "Fura2Ca")
    fura2 = np.array(fura2)
    
    mean = fura2[:, :int(3000/dt)].mean(axis=2)
    fura_scaled = (fura2 - mean[:, :, None])/mean[:, :, None]
    fluo_vals = fura_scaled.max(axis=2)
    return fluo_vals, voxels

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
RyR_dir = "Ca_wave_simple_SERCA_no_SOCE_Breit_2018"
output_name = "all"
stims = ["0350", "0700", "1050"]
stim_labels = ["1 uM Ca injection", "2.5 uM Ca injection", "5 uM Ca injection"]
branch_diams = [1.2, 2.4, 6.0]
Fura_specie = "Fura2Ca"
t_stim = 3000  # sec
dye_base = "model_RyR2CaM_simple_SERCA_SOCE_Fura2_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
dye_base_RyR = "model_RyR_simple_SERCA_Fura2_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
nodye_base = "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(1, 3, figsize=(15, 10))

#  Dye figs

for i, x in enumerate(ax):
    x.set_title(stim_labels[i], fontsize=20)
    for k, b_diam in enumerate(branch_diams):
        fname = dye_base % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname)
        print(full_name)
        try:
            fluo_vals, voxels = get_fluo(full_name)
            max_fluo_vals = fluo_vals.mean(axis=0)
            max_fluo_error = fluo_vals.std(axis=0)/fluo_vals.shape[0]**0.5
        
            vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
            x.plot(vox_axis[25:-25], max_fluo_vals[25:-25], color=colors_ctrl[k], 
                   label="dend diam=%2.1f um" % b_diam)
            x.fill_between(vox_axis[25:-25], max_fluo_vals[25:-25]-max_fluo_error[25:-25],
                           max_fluo_vals[25:-25]+max_fluo_error[25:-25],
                           color=colors_dye[k], alpha=0.2)
        except ValueError:
            pass
        
        

        fname = dye_base_RyR % (b_diam, stims[i])
        full_name = os.path.join(cur_dir, RyR_dir, fname)
        try:
            fluo_vals, voxels = get_fluo(full_name)
            max_fluo_vals = fluo_vals.mean(axis=0)
            max_fluo_error = fluo_vals.std(axis=0)/fluo_vals.shape[0]**0.5
            
            vox_axis = np.linspace(-voxels[-1]/2, voxels[-1]/2, len(voxels))
            x.plot(vox_axis[25:-25], max_fluo_vals[25:-25], color=colors_dye[k], 
                   label="dend diam=%2.1f um no CaM" % b_diam)
            x.fill_between(vox_axis[25:-25], max_fluo_vals[25:-25]-max_fluo_error[25:-25],
                           max_fluo_vals[25:-25]+max_fluo_error[25:-25],
                           color=colors_dye[k], alpha=0.2)
        except TypeError:
            pass

        
        if not i:
            x.legend()
            x.set_ylabel("%Fluorescence", fontsize=20)
        else:
            x.set_yticklabels([])
    
        x.tick_params(labelsize=14)
        x.set_xlabel("Distance from stim site [um]", fontsize=20)


mini, maxi = get_ylim(ax)
set_ylim(ax, mini, maxi)
fig.savefig("Ca_dye_effects.eps", dpi=100, bbox_inches="tight")
fig.savefig("Ca_dye_effects.png", dpi=100, bbox_inches="tight")


