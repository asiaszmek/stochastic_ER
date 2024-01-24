import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {"tubes_3s_injection": "tab:blue",
          "tubes": "tab:green",
          "baloon_3s_injection": "tab:cyan",
          "baloon": "tab:olive",
}
labels = {"no_SOCEtubes": "no SOCE RyR dis-inh.",
          "tubes": "RyR dis-inh.",
          "no_SOCEbaloon": "no SOCE",
          "baloon": "ctrl",
}

def set_ylim(ax, mini, maxi):
    for x in ax:
        x.set_ylim([mini + 0.01, maxi + 0.05])


reg_list = ["dend25", "dend26", "dend27"]

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
t_start = 3000
idx_start = t_start

base_noSOCE = "model_RyR%s_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
base_SOCE = "model_RyR%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

CaM_noSOCE = "model_RyR2CaM%s_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
fname = {
    "baloon": CaM_SOCE,
    "tubes": base_SOCE
    }
directory ={
    "baloon": "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "tubes": "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    } 
stim_types = ["_3s_injection", ""]
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))
suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}



for i, x in enumerate(ax):
    mini = []
    maxi = []
    for k, b_diam in enumerate(branch_diams):
        x[k].set_title(stim_labels[i]+ " diam %2.1f um" % b_diam,
                       fontsize=12)

        for stim in stim_types:
            for inh in ["baloon", "tubes"]:
                fname_SOCE = fname[inh] % (stim, b_diam, stims[i])
                full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
                try:
                    voxels, time, ca = utils.get_conc(full_name, ["Ca"],
                                                      reg_list,
                                                      output_name)
                except FileNotFoundError and OSError:
                    print("Could not find", full_name)
                    continue
                output = ca.mean(axis=1)
                x[k].plot(time/1000, output, colors[inh+stim],
                          label=labels[inh]+" "+suffix[stim])
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
            x[k].set_xlabel("Time [sec]", fontsize=20)
            x[k].tick_params(labelsize=14)

    set_ylim(x, min(mini), max(maxi))
ax[2, 2].legend()

fig.savefig("Ca_decay_stim_point.eps", dpi=100, bbox_inches="tight")
plt.show()
