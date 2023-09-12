import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import sys
from scipy.constants import Avogadro
import utility_functions as utils

NA = Avogadro*1e-23


if __name__ == '__main__':
    specie_list = ["Ca"]
    specie = "Ca"
    reg_list = ["dend", "dend01", "dend02", "dend03", "dend04", "dend05", "dend06",
                 "dend07", "dend08", "dend09", "dend10", "dend11"]
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    for fname in sys.argv[1:]:
        my_file = h5py.File(fname, 'r')
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
        vmax = 0

        for key in conc_dict:
            new_max = conc_dict[key].max()
            if new_max > vmax:
                vmax = new_max
        for key in conc_dict:
            fig, ax = plt.subplots(1, 1)
            time = time_dict[key]
            im = ax.imshow(conc_dict[key].T, aspect="auto",
                           interpolation="none",
                           origin="lower", extent = [time[0]*1e-3,
                                                     time[-1]*1e-3,
                                                     voxels[0],
                                                     voxels[-1]],
                           cmap=plt.get_cmap("Reds"), vmin=0, vmax=vmax)
            fig.colorbar(im)
            ax.set_title("%s %s" % (specie, trial))
    plt.show()
                          
