import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils



colors = ["tab:blue", "tab:olive", "tab:green"]

def set_ylim(ax, mini, maxi):
    for x in ax:
        x.set_ylim([mini + 0.01, maxi + 0.05])


reg_list = ["dend25", "dend26", "dend27"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
ctrl_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
aging_dir = "Ca_wave_aging"
output_name = "all"
stim = "0350"
stim_label = "2 uM Ca injection"
branch_diams = [1.2, 2.4, 6.0]
t_start = 3000
idx_start = t_start
labels = {"aging": "old age",
          "ctrl": "ctrl",
}

base_ctrl = "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

base_aging = "model_aging%s_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fname = {
    "ctrl": base_ctrl,
    "aging": base_aging
    }
directory ={
    "ctrl": ctrl_dir,
    "aging": aging_dir,
    } 
stim_types = ["_3s_injection", ""]

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))


for i, b_diam in enumerate(branch_diams):
    for inh in ["ctrl", "aging"]:
        for k, s in enumerate(stim_types):
            ax[k].set_title(stim_label+ " diam %2.1f um" % b_diam,
                            fontsize=12)
            fname_SOCE = fname[inh] % (s, b_diam, stim)
            full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
            try:
                voxels, time, ca = utils.get_conc(full_name, ["Ca"],
                                                      reg_list,
                                                      output_name)
            except FileNotFoundError and OSError:
                print("Could not find", full_name)
                continue
            output = ca.mean(axis=1)
            max_idx = output.argmax()
            if inh == "aging":
                ax[k].plot(time[max_idx:len(time)//2:100]/1000,
                           output[max_idx:len(time)//2:100], color=colors[i],
                           label=labels[inh]+" diam "+ str(b_diam),
                           linestyle="", marker="d",  fillstyle="none")
            else:
                ax[k].plot(time[max_idx:len(time)//2:100]/1000,
                           output[max_idx:len(time)//2:100], color=colors[i],
                           label=labels[inh]+" diam "+ str(b_diam),
                           linestyle="", marker="d",  fillstyle="full")
            
ax[0].set_ylabel("Calcium [uM]", fontsize=20)
ax[0].set_xlabel("Time [sec]", fontsize=20)
ax[1].set_xlabel("Time [sec]", fontsize=20)
ax[0].tick_params(labelsize=14)
ax[1].tick_params(labelsize=14)
ax[1].legend()
ax[0].set_title("3 ms stim. duration")
ax[1].set_title("40 ms stim. duration")
ylim = max(ax[0].get_ylim() + ax[1].get_ylim())
ax[0].set_ylim([0, ylim])
ax[1].set_ylim([0, ylim])
fig.savefig("Temporal_ctrl_aging_350.eps", dpi=100, bbox_inches="tight")
fig.savefig("Temporal_ctrl_aging_350.png", dpi=100, bbox_inches="tight")
plt.show()
