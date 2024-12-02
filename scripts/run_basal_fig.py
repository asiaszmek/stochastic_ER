import os
import sys
import h5py
import numpy as np
from scipy.fft import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
import utility_functions as ut
colors =  {
    "1.2": 'tab:blue',
    "2.4": 'tab:purple',
    "6.0": 'tab:green'
}

names_dict = {
    "normal PMCA + RyR2CaM" : "model_RyRCaM_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    "normal PMCA + no RyR2": "model_noRyR_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    "normal PMCA + RyR2": "model_RyR_simple_SERCA_SOCE_tubes_diam_%s_um_2_um_dendrite.h5",
    "low PMCA + RyR2CaM": "model_RyRCaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "low PMCA, RyR2 no CaM": "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "low PMCA + 50% RyR2 + 50% RyR2CaM": "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5"
}

dend_diam = ["1.2", "2.4", "6.0"]
trial = "trial0"
output = "__main__"
if __name__ == "__main__":
    data_dir = os.path.join("..", "stacked_ER")
    x_labels = list(names_dict.keys())
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))
    for i, d in enumerate(dend_diam):
        means = []
        stds = []
        for fname in names_dict.values():
            path = os.path.join(data_dir, fname % d)
            print(path)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            data = ut.get_populations(my_file, trial=trial, output=output)
            specie_idx = ut.get_all_species(my_file, output=output).index("Ca")
            volume = sum(list(ut.region_volumes(my_file).values()))

            conc = ut.nano_molarity(data[:, :, specie_idx].sum(axis=1), volume)
            means.append(conc[len(conc)//2:].mean())
            stds.append(conc[len(conc)//2:].std())
        ax1[i].errorbar(y=x_labels, x=means, xerr=stds,
                        color=colors[d],
                        marker="o", linestyle="")
        ax1[i].set_title("diam %s um" % d)
        ax1[i].set_yticklabels([])
    ax1[0].set_yticklabels(x_labels)
    mini = min([min(x.get_xlim()) for x in ax1])
    maxi = max([max(x.get_xlim()) for x in ax1])
    for x in ax1:
        x.set_xlabel("Basal Ca [nM]")
        x.set_xlim([mini, maxi])
    fig1.savefig("mean_basal_ca.png", dpi=100,
                 bbox_inches="tight")
                    

