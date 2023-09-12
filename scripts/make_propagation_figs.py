import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import utility_functions as utils

if __name__ == '__main__':
    reg_list = ["dend", "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09", "dend10",
                "dend11", "dend12", "dend13", "dend14", "dend15",
                "dend16", "dend17", "dend18", "dend19", "dend20",
                "dend21", "dend22", "dend23", "dend24", "dend25",
                "dend26", "dend27", "dend28", "dend29", "dend30",
                "dend31", "dend32", "dend33", "dend34", "dend35",
                "dend36", "dend37", "dend38", "dend39", "dend40",
                "dend41", "dend42", "dend43", "dend44", "dend45",
                "dend46", "dend47", "dend48", "dend49", "dend50",
                "dend51"]
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    for fname in sys.argv[1:]:
        my_file = h5py.File(fname, 'r')
        specie_list = utils.get_all_species(my_file, output="all")
        conc_dict = {}
        times = []
        voxels = {}
        for trial in my_file.keys():
            if trial == "model":
                continue
            for specie in specie_list:
                conc, voxels[specie] = utils.get_dynamics_in_region(my_file,
                                                                    specie,
                                                                    reg_list, trial,
                                                                    "all")
                
                if specie not in conc_dict:
                    conc_dict[specie] = list()
                    conc_dict[specie].append(conc)
            times.append(get_times(my_file, trial, "all"))
        for specie in conc_dict.keys():
            amin = np.argmin([len(c) for c in conc_dict[specie]])
            shape = conc_dict[specie][amin].shape
            avg_conc = np.zeros(shape)
            for c in conc_dict[specie]:
                avg_conc += c[:shape[0], :]
            avg_conc = avg_conc/len(conc_dict[specie])
            length = voxels[specie][-1]-voxels[specie][0]
            len_array = np.linspace(voxels[specie][0], voxels[specie][-1],
                                    shape[1]//2)
            
            fig, ax = plt.subplots(1, 1)
            for i, x in enumerate(len_array):
                ax.plot(times[amin], avg_conc[:, 2*i:2*i+2].mean(axis=1),
                        marker=".", label="%4.2f um" % x, linestyle='None')
            ax.set_title("%s %s" % (fname, specie))
            if specie == "Ca":
                ax.set_yscale("log")
            ax.legend()
    plt.show()
                          
