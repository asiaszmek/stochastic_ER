import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = ["tab:blue", "tab:olive", "tab:green"]

labels = {"tubes": "uniformly distributed RyR2CaM",
          "baloon": "RyR2CaM in the dendritic membrane",
}


reg_list = ["dend25", "dend26", "dend27"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
basic_RyR_SOCE_dir = "Ca_wave_simple_SERCA_SOCE"
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stim = "0350"
stim_label = "2 uM Ca injection"
branch_diams = [1.2, 2.4, 6.0]

t_start = 3000
idx_start = t_start

base_SOCE = "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM%s_simple_SERCA_SOCE_baloon_diam_%2.1f_um_50_um_%s_nM.h5"

fname = {
    "baloon": CaM_SOCE,
    "tubes": base_SOCE
    }

directory ={
    "baloon": "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "tubes": "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    }

stim_types = ["_3s_injection", ""]
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}


for i, b_diam in enumerate(branch_diams):
    for inh in ["tubes", "baloon"]:
        for k, s in enumerate(stim_types):
            fname_SOCE = fname[inh] % (s, b_diam, stim)
            full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
            try:
                voxels, time, ca = utils.get_conc(full_name, ["Ca"],
                                                  reg_list,
                                                  output_name)
            except TypeError:
                print("Could not find", full_name)
                continue
            output = ca.mean(axis=1)
            max_idx = output.argmax()
            if inh == "tubes":
                ax[k].plot(time[max_idx:len(time)//2:50]/1000,
                           output[max_idx:len(time)//2:50], color=colors[i],
                           label=labels[inh]+" diam "+ str(b_diam),
                           linestyle="", marker="d",  fillstyle="none")
            else:
                ax[k].plot(time[max_idx:len(time)//2:50]/1000,
                           output[max_idx:len(time)//2:50], color=colors[i],
                           label=labels[inh]+" diam "+ str(b_diam),
                           linestyle="", marker="d",  fillstyle="full")

ax[0].set_ylabel("Calcium [uM]", fontsize=20)
ax[0].set_xlabel("Time [sec]", fontsize=20)
ax[1].set_xlabel("Time [sec]", fontsize=20)
ax[0].tick_params(labelsize=14)
ax[1].tick_params(labelsize=14)
ax[1].legend(loc='upper left', bbox_to_anchor=(1, 0.5))
ax[0].set_title("3 ms stim. duration")
ax[1].set_title("40 ms stim. duration")
ylim = max(ax[0].get_ylim() + ax[1].get_ylim())
ax[0].set_ylim([0, ylim])
ax[1].set_ylim([0, ylim])

fig.savefig("Ca_decay_stim_3_ms_40_ms_temp_comp.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Ca_decay_stim_3_ms_40_ms_temp_comp.png", dpi=100,
            bbox_inches="tight")

plt.show()
