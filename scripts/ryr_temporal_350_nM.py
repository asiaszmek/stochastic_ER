import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {"_no_CaM_3s_injection": "tab:blue",
          "_no_CaM": "tab:green",
          "_CaM_3s_injection": "tab:cyan",
          "_CaM": "tab:olive",
}
labels = {"no_SOCE_no_CaM": "no SOCE RyR dis-inh.",
          "_no_CaM": "RyR dis-inh.",
          "no_SOCE_CaM": "no SOCE",
          "_CaM": "ctrl",
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

stim_types = ["_3s_injection", ""]
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}


for k, b_diam in enumerate(branch_diams):
    for inh in ["_no_CaM", "_CaM"]:
        for k, s in enumerate(stim_types):
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
            print(max(output))
            ax[k].plot(time/1000, output, 
                    label=labels[inh]+" diam "+ str(b_diam))

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

fig.savefig("Ca_decay_stim_3_ms_40_ms_temp_comp.svg", dpi=100, bbox_inches="tight")
plt.show()
