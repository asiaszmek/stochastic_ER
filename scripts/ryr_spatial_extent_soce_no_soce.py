import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {1.2: "tab:blue",
          2.4: "tab:green",
          6.0: "tab:olive",
}
labels = {"_SOCE": "control",
          "_no_SOCE": "control no SOCE",
}


reg_list = ["dend25", "dend26", "dend27"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
stims = ["0175", "0350", "0700", "1050", "2000"]
branch_diams = [1.2, 2.4, 6.0]


t_start = 3000
idx_start = t_start

base_SOCE = "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM%s_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fname = {
    "_SOCE": CaM_SOCE,
    "_no_SOCE": base_SOCE
    }

directory ={
    "_SOCE": "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "_no_SOCE": "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
    }

stim_types = ["_3s_injection", ""]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}


for k, b_diam in enumerate(branch_diams):
    for inh in ["_SOCE", "_no_SOCE"]:
        for j, s in enumerate(stim_types):
            y = []
            x = []
            for stim in stims:
                fname_SOCE = fname[inh] % (s, b_diam, stim)
                full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
                try:
                    voxels, time, ca = utils.get_conc(full_name, ["Ca"],
                                                      reg_list,
                                                      output_name)
                except TypeError:
                    print("Could not find", full_name)
                    continue
                y.append(np.max(ca.mean(axis=1))/1000)
                x.append(injections[b_diam][stim][s]/b_diam)
            if inh == "_SOCE" and s == "":
                ax.plot(x, y,color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" 40 ms",
                        linestyle="")
            elif inh == "_no_SOCE" and s == "":
                ax.plot(x, y, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" 40 ms",
                        linestyle="",
                        fillstyle="none")
            elif inh == "_SOCE" and s == "_3s_injection":
                ax.plot(x, y, color=colors[b_diam], marker="o",
                        label=labels[inh]+" diam "+str(b_diam)+" 3 ms",
                        linestyle="")
            else:
                ax.plot(x, y, color=colors[b_diam], marker="o",
                        label=labels[inh]+" diam "+str(b_diam)+" 3 ms",
                        linestyle="", fillstyle="none")

ax.set_ylabel("max(Calcium) [uM]", fontsize=20)
ax.set_xlabel("total injected ions/diam [1/um]", fontsize=20)
ax.tick_params(labelsize=14)
ax.legend(loc='lower left', bbox_to_anchor=(1, 0.5))


fig.savefig("Spatial_sum_soce.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Spatial_sum_soce.png", dpi=100,
            bbox_inches="tight")

plt.show()
