import os
import utility_functions as utils


colors = {1.2: "tab:blue",
          2.4: "tab:purple",
          6.0: "tab:green",
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

base_SOCE = "model_RyR%s_simple_SERCA_SOCE_baloon_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM%s_simple_SERCA_SOCE_baloon_diam_%2.1f_um_50_um_%s_nM.h5"

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

suffix = {"": " EPSP",
          "_3s_injection": "bAP"}


for k, b_diam in enumerate(branch_diams):
    for inh in ["_no_CaM", "_CaM"]:
        for j, s in enumerate(stim_types):
            y = []
            y_err = []
            x = []
            for stim in stims:
                fname_SOCE = fname[inh] % (s, b_diam, stim)
                full_name = os.path.join(cur_dir, directory[inh],
                                         fname_SOCE)
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
            if inh == "_CaM" and s == "":
                ax.errorbar(x, y, yerr=y_err,color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)
                        +" EPSP",
                        linestyle="")
            elif inh == "_no_CaM" and s == "":
                ax.errorbar(x, y, yerr=y_err, color=colors[b_diam], marker="d",
                        label=labels[inh]+" diam "+str(b_diam)+" EPSP",
                        linestyle="",
                        fillstyle="none")
            elif inh == "_CaM" and s == "_3s_injection":
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
ax.legend()


fig.savefig("Spatial_sum_baloon.eps", dpi=100,
            bbox_inches="tight")
fig.savefig("Spatial_sum_tubes.png", dpi=100,
            bbox_inches="tight")


