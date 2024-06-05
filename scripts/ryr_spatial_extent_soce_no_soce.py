import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {1.2: "tab:blue",
          2.4: "tab:purple",
          6.0: "tab:green",
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

suffix = {"": " EPSP",
          "_3s_injection": " bAP"}


for k, b_diam in enumerate(branch_diams):
    for inh in ["_SOCE", "_no_SOCE"]:
        for j, s in enumerate(stim_types):
            y = []
            x = []
            y_err = []
            for stim in stims:
                fname_SOCE = fname[inh] % (s, b_diam, stim)
                full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
                ca_dict, time_dict = utils.get_conc(full_name, ["Ca"],
                                                    reg_list,
                                                    output_name)
                ca_out = [ca_dict["Ca"][key]
                          for key in ca_dict["Ca"].keys()]
                ca = np.array(ca_out)/1000
            
                output = ca.max(axis=0)
                y.append(np.mean(output))
                y_err.append(np.std(output)/len(output)**0.5)
                x.append(injections[b_diam][stim][s]/b_diam)
            if inh == "_SOCE" and s == "":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" EPSP",
                        linestyle="")
            elif inh == "_no_SOCE" and s == "":
                ax.errorbar(x, y, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" EPSP",
                        linestyle="",
                        fillstyle="none")
            elif inh == "_SOCE" and s == "_3s_injection":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="o",
                        label=labels[inh]+" diam "+str(b_diam)+" bAP",
                        linestyle="")
            else:
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="o",
                        label=labels[inh]+" diam "+str(b_diam)+" bAP",
                        linestyle="", fillstyle="none")

ax.set_ylabel("max(Calcium) [uM]", fontsize=20)
ax.set_xlabel("total injected ions/diam [1/um]", fontsize=20)
ax.tick_params(labelsize=14)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


fig.savefig("Spatial_sum_soce.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Spatial_sum_soce.png", dpi=100,
            bbox_inches="tight")

plt.show()
