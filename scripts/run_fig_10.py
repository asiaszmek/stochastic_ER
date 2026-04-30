import os
import sys
import h5py
from matplotlib.patches import Rectangle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import utility_functions as ut
plt.rcParams['text.usetex'] = True
colors =  {
    "1.2": 'tab:blue',
    "2.4": 'tab:purple',
    "6.0": 'tab:green'
}

names_dict = {
   
    "80\%\n100\%\n0\%":
    os.path.join("model_RyRCaM_0.8_PMCA",
                 "model_RyRCaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5"),
    "80\%\n50\%\n50\%":
    os.path.join("model_RyR_RyRCaM_0.8_PMCA",
                 "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5"),
    "80\%\n0\n100\%":
    os.path.join("model_RyR_0.8_PMCA",
                 "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5"),
 
    "80\%\n100\%\n100\%":
    os.path.join("model_2x_RyR_RyRCaM_0.8_PMCA",
                  "model_2x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5"),
    "80\%\n0\n200\%":
    os.path.join("model_2xRyR_0.8_PMCA",
                 "model_2xRyR_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_2_um_dendrite.h5"),
    
    
}

dend_diam = ["1.2", "2.4", "6.0"]
output = "__main__"

        

def get_ca_conc(my_file, trial, output=output):
    try:
        data = ut.get_populations(my_file, trial=trial,
                                  output=output)
    except KeyError:
        return
                
    specie_idx = ut.get_all_species(my_file,
                                    output=output).index("Ca")

    volume = np.array([g[12] for g in grid_list])
    conc =  ut.nano_molarity(data[:, :, specie_idx],
                                         volume)
    new_conc = np.reshape(conc, (data.shape[0],
                                 102,
                                 len(volume)//102)).mean(axis=(1,2))
              
    return new_conc


dt = 0.2 # s

if __name__ == "__main__":
    data_dir = ".."
    x_labels = list(names_dict.keys())
    fig_wave_f, ax_wave_f = plt.subplots(1, 1,
                                         figsize=(1*7, 5))
    for i, d in enumerate(dend_diam):
        f_mean = []
        f_std = []
        xlabels = []
        for key, fname in names_dict.items():
            path = os.path.join(data_dir, fname % d)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            peak_no = []
            length = []
            for trial in ["trial0", "trial1", "trial2", "trial3", "trial4"]:
                new_conc = get_ca_conc(my_file, trial)

                if new_conc is None:
                    continue
                peak = 0
                where = np.where(new_conc > 2.5*76)
                if len(where[0]):
                    prev = where[0][0]
                    peak = 1
                    for idx in where[0][1:]:
                        if idx >= prev +1 and idx <= prev + 20:
                            prev = idx
                        else:
                            peak += 1
                            prev = idx
                peak_no.append(peak)
                length.append(len(new_conc)*dt)
            xlabels.append(key)
            m_f = (np.array(peak_no)/np.array(length)).mean()
            std_f = (np.array(peak_no)/np.array(length)).std()/(len(peak_no))**0.5
            f_mean.append(m_f)
            f_std.append(std_f)

        ax_wave_f.errorbar(x=xlabels, y=f_mean,
                            yerr=f_std,
                            color=colors[d],
                            marker="o", linestyle="",
                            label=r"%s $\mathrm{\mu m}$ diam" % d)

                                
    ax_wave_f.set_ylabel(r"$\mathrm{Ca_i}$ wave frequency (Hz)",
                          fontsize=15)
        
    ax_wave_f.set_xticklabels(xlabels)
    ax_wave_f.tick_params(axis='x', labelsize=15)
    ax_wave_f.tick_params(axis='y', labelsize=15)
    legend = "PMCA kcat\nRyR2CaM\nRyR2"#\nSOCE"
    ax_wave_f.legend()
    rect = Rectangle((2.5,-0.001), 1, 0.01, edgecolor="r", facecolor="none")
    ax_wave_f.add_patch(rect)
    ax_wave_f.text(2.5,0.009901, "AD model", color="r")
    
    ax_wave_f.text(-1.25, min(ax_wave_f.get_ylim())
                    -(max(ax_wave_f.get_ylim())
                      -min(ax_wave_f.get_ylim()))*0.1698,
                           legend,
                    horizontalalignment='left', fontsize=15)
    

    fig_wave_f.savefig("mean_wave_frequency.png", dpi=100,
                       bbox_inches="tight")
    fig_wave_f.savefig("mean_wave_frequency.eps", dpi=100,
                       bbox_inches="tight")
