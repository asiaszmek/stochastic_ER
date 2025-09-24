import os
import h5py
import utility_functions as utils
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True


if __name__ == '__main__':
    specie_list = ["Ca"]
    specie = "Ca"
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
    fnames = [
        os.path.join("..",
                     "model_RyR_RyRCaM_0.8_PMCA",
                     "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_2.4_um_2_um_dendrite.h5"),
        os.path.join("..",
                     "model_RyR_0.8_PMCA",
                     "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_2.4_um_10_um_dendrite.h5")
    ]
    figs, axes = [], []
    for fname in fnames:
        my_file = h5py.File(fname, 'r')
        new_fname = os.path.split(fname)[-1]
        conc_dict = {}
        time_dict = {}
        for trial in my_file.keys():
            if trial == "model":
                continue
            conc, voxels = utils.get_dynamics_in_region(my_file,
                                                        specie_list,
                                                        reg_list, trial, "__main__")
            conc_dict[trial] = conc
            time = utils.get_times(my_file, trial, "__main__")
            time_dict[trial] = time
        vmin = 0
        vmax = 1200
        diam = fname.split("diam_")[-1][:3]
        for key in conc_dict:
            fig, ax = plt.subplots(1, 1)
            time = time_dict[key]
            im = ax.imshow(conc_dict[key].T, aspect="auto",
                           interpolation="none",
                           origin="lower", extent = [time[0]*1e-3,
                                                     time[-1]*1e-3,
                                                     voxels[0],
                                                     voxels[-1]],
                           cmap=plt.get_cmap("Reds"))
            ax.set_xlabel(r"time (s)", fontsize=15)
            ax.set_ylabel(r"dendrite $(\mathrm{\mu m})$", fontsize=15)
            fig.colorbar(im)
            
            ax.set_title(r"%s dynamics in %s $\mathrm{\mu m}$ dend" % (specie, diam),
                         fontsize=14)
            fig.savefig(new_fname[:-3]+"_"+key+".png", dpi=100,
                        bbox_inches="tight")
            fig.savefig(new_fname[:-3]+"_"+key+".eps", dpi=100,
                        bbox_inches="tight")
    
    
    
 
                          
