import os
import utility_functions as utils
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

legend_elements = [Line2D([0], [0], color='k', marker="o", fillstyle="full",
                          lw=0, label='ctrl'),
                   Line2D([0], [0], color="k", marker='o', fillstyle="none",
                          lw=0, label="no SOCE"),
                  ]

colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}


directories = [
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_0700_nM.h5"
    ],
    [
        "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
        "model_RyR2CaM_simple_SERCA_tubes_diam_2.4_um_50_um_0700_nM.h5"
    ],
]


types = [ "RyR2CaM+SOCE", "RyR2CaM"]
marker = ["o", "o"]
fillstyle = ["full", "none"]
stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = [ "2.4"]
base = "dend"
reg_list = [base, "dend01", "dend02", "dend03", "dend04",
            "dend05", "dend06", "dend07", "dend08", "dend09",]
for i in range(10, 102, 1):
    reg_list.append("%s%d" %(base, i))
linestyles = ["solid", "dotted"]
trial = "trial1"
if __name__ == '__main__':
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    for k, (d, fname) in enumerate(directories):
        my_path = os.path.join("..", d, fname)
        conc_dict, times_dict = utils.get_conc(my_path,
                                               ["CaER"],
                                               reg_list,
                                               "all")
        length = len(times_dict[trial])
        print(conc_dict["CaER"][trial].min())
        ca_er = conc_dict["CaER"][trial].sum(axis=0)

        ax1.plot(times_dict[trial]/1000, ca_er, label=types[k], color="tab:purple", linestyle=linestyles[k])
        print(min(ca_er))
    ax1.set_xlabel("time (s)", fontsize=15)
    ax1.set_ylabel("CaER $(\mu\mathrm{M})$", fontsize=15)
    ax1.tick_params(axis='x', labelsize=15)
    ax1.tick_params(axis='y', labelsize=15)
    ax1.legend()
    fig1.savefig("min_CaER_SOCE_no_SOCE.png", dpi=100,
                 bbox_inches="tight")
    fig1.savefig("min_CaER_SOCE_no_SOCE.eps", dpi=100,
                bbox_inches="tight")
