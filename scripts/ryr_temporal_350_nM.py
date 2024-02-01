import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = ["tab:blue", "tab:purple", "tab:green"]
colors_aging = ["tab:cyan", "tab:pink", "tab:olive"]
labels = {"_no_CaM": "no CaM",
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
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}


for i, b_diam in enumerate(branch_diams):
    for inh in ["_no_CaM", "_CaM"]:
        for k, s in enumerate(stim_types):
            fname_SOCE = fname[inh] % (s, b_diam, stim)
            full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
            ca_dict, time_dict = utils.get_conc(full_name, ["Ca"],
                                                reg_list,
                                                output_name)
    
            ca_out = [ca_dict["Ca"][key].mean(axis=0) for key in ca_dict["Ca"].keys()]
            try:
                time = time_dict["trial0"]
            except KeyError:
                continue
            
            ca = np.array(ca_out)/1000
            output = ca.mean(axis=0)
            ca_err = ca.std(axis=0)/10**0.5
            max_idx = output.argmax()
            print(max(ca_err))
            if inh == "_no_CaM":
                ax[k].plot(time[max_idx:len(time)//4]/1000,
                           output[max_idx:len(time)//4], color=colors_aging[i],
                           label=labels[inh]+" diam "+ str(b_diam))
                ax[k].fill_between(time[max_idx:len(time)//4]/1000,
                                   output[max_idx:len(time)//4]
                                   -ca_err[max_idx:len(time)//4],
                                   output[max_idx:len(time)//4]
                                   +ca_err[max_idx:len(time)//4],
                                   color=colors_aging[i],
                                   alpha=0.2)
            else:
                ax[k].plot(time[max_idx:len(time)//4]/1000,
                           output[max_idx:len(time)//4], color=colors[i],
                           label=labels[inh]+" diam "+ str(b_diam))
                ax[k].fill_between(time[max_idx:len(time)//4]/1000,
                                   output[max_idx:len(time)//4]
                                   -ca_err[max_idx:len(time)//4],
                                   output[max_idx:len(time)//4]
                                   +ca_err[max_idx:len(time)//4],
                                   color=colors[i], alpha=0.2)

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

fig.savefig("Ca_decay_stim_3_ms_40_ms_temp_comp_CaM.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Ca_decay_stim_3_ms_40_ms_temp_comp_CaM.png", dpi=100,
            bbox_inches="tight")

plt.show()
