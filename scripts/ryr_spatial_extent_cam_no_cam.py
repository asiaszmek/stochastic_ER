import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {1.2: "tab:blue",
          2.4: "tab:green",
          6.0: "tab:purple",
}
labels = {"no_SOCE_no_CaM": "no CaM no SOCE",
          "_no_CaM": "no CaM",
          "no_SOCE_CaM": "ctrl no SOCE",
          "_CaM": "ctrl",
}


reg_list = ["dend25", "dend26", "dend27"]

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
basic_RyR_SOCE_dir = "Ca_wave_simple_SERCA_SOCE"
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stims = ["0175", "0350", "0700", "1050", "2000"]
stim_label = "2 uM Ca injection"
branch_diams = [1.2, 2.4, 6.0]

injections = {
    1.2:{
        "0175":{
            "":20*1200,
            "_3s_injection":1.5*12000,
        },

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
         "2000":{
            "":40*4800,
            "_3s_injection":3*48000,
        },
    },
    2.4:{
        "0175":{
            "":20*2000,
            "_3s_injection":1.5*20000,
        },

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
        "2000":{
            "":40*8000,
            "_3s_injection":3*80000,
        },
    },
    6.0:{
        "0175":{
            "":20*4000,
            "_3s_injection":1.5*40000,
        },

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
        "2000":{
            "":40*16000,
            "_3s_injection":3*160000,
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
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

suffix = {"": "40 ms stim",
          "_3s_injection": "3 ms stim"}


for k, b_diam in enumerate(branch_diams):
    for inh in ["_no_CaM", "_CaM"]:
        for j, s in enumerate(stim_types):
            y = []
            y_err = []
            x = []
            for stim in stims:
                fname_SOCE = fname[inh] % (s, b_diam, stim)
                full_name = os.path.join(cur_dir, directory[inh], fname_SOCE)
                ca_dict, time_dict = utils.get_conc(full_name, ["Ca"],
                                                    reg_list,
                                                    output_name)
                try:
                    ca_out = utils.get_array(ca_dict, "Ca")
                except ValueError:
                    continue
                ca = np.array(ca_out)/1000
                ca = ca.mean(axis=1)
                output = ca.max(axis=1)
                y.append(np.mean(output))
                x.append(injections[b_diam][stim][s]/b_diam)
                y_err.append(output.std()/len(output)**0.5)
            print(x, y, y_err)
            if inh == "_CaM" and s == "":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" 40 ms",
                        linestyle="")
            elif inh == "_no_CaM" and s == "":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" 40 ms",
                        linestyle="",
                        fillstyle="none")
            elif inh == "_CaM" and s == "_3s_injection":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="o",
                        label=labels[inh]+" diam "+str(b_diam)+" 3 ms",
                        linestyle="")
            else:
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="o",
                        label=labels[inh]+" diam "+str(b_diam)+" 3 ms",
                        linestyle="", fillstyle="none")

ax.set_ylabel("max(Calcium) [uM]", fontsize=20)
ax.set_xlabel("total injected ions/diam [1/um]", fontsize=20)

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.tick_params(axis='both', which='major', labelsize=14)

fig.savefig("Spatial_sum_cam.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Spatial_sum_cam.png", dpi=100,
            bbox_inches="tight")

plt.show()
