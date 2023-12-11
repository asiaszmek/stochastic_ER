import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {1.2: "tab:blue",
          2.4: "tab:green",
          6.0: "tab:olive",
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
stims = ["0350", "0700", "1050"]
stim_label = "2 uM Ca injection"
branch_diams = [1.2, 2.4, 6.0]

injections = {
    1.2:{
        "0350":{
            "":40*1200,
            "_3s_injection":3*12000,
        },
        "0700":{
            "":40*2400,
            "_3s_injection":3*24000,
        },
        "1050":{
            "":40*4000,
            "_3s_injection":3*40000,
        },
    },
    2.4:{
        "0350":{
            "":40*2000,
            "_3s_injection":3*20000,
        },
        "0700":{
            "":40*4000,
            "_3s_injection":3*40000,
        },
        "1050":{
            "":40*6000,
            "_3s_injection":3*60000,
        },

    },
    6.0:{
        "0350":{
            "":40*4000,
            "_3s_injection":3*40000,
        },
        "0700":{
            "": 40*8000,
            "_3s_injection":3*80000,
        },
        "1050":{
            "":40*12000,
            "_3s_injection":3*120000,
        },

    },
}

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
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}


for k, b_diam in enumerate(branch_diams):
    for inh in ["_no_CaM", "_CaM"]:
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
                except FileNotFoundError and OSError:
                    print("Could not find", full_name)
                    continue
                y.append(np.max(ca.mean(axis=1))/1000)
                x.append(injections[b_diam][stim][s]/b_diam)
            if inh == "_CaM" and s == "":
                ax.plot(x, y,color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" 40 ms",
                        linestyle="")
            elif inh == "_no_CaM" and s == "":
                ax.plot(x, y, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" 40 ms",
                        linestyle="",
                        fillstyle="none")
            elif inh == "_CaM" and s == "_3s_injection":
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
ax.legend()


fig.savefig("Spatial_sum.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Spatial_sum.png", dpi=100,
            bbox_inches="tight")

plt.show()
