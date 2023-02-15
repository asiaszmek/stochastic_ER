import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import sys
from scipy.constants import Avogadro

NA = Avogadro*1e-23


def nano_molarity(N, V):
    return 10 * N / V / NA


def get_grid_list(My_file):
    return np.array(My_file['model']['grid'])


def get_times(My_file, trial='trial0', output="__main__"):
    return np.array(My_file[trial]['output'][output]['times'])


def get_outputs(my_file):
    return my_file['model']['output'].keys()


def get_populations(my_file, trial='trial0', output='__main__'):
    return np.array(my_file[trial]['output'][output]['population'])


def get_all_species(My_file, output="__main__"):
    return [s.decode('utf-8') for s in My_file['model']['output'][output]['species']]


def get_dend_indices(grid, region="dend"):
    out = {}
    volumes = {}
    if not isinstance(region, list):
        region = [region]
    for i, line in enumerate(grid):
        if line[15].decode('utf-8') in region:
            pos = abs(np.round(line[0], 3))
            if pos in out:
                out[pos].append(i)
            else:
                out[pos] = [i]
            if pos in volumes:
                volumes[pos] += line[12]
            else:
                volumes[pos] = line[12]
    return out, volumes


def get_dynamics_in_region(my_file, specie, region, trial,
                           output="__main__"):
    if not isinstance(specie, list):
        specie = [specie]
    my_grid = get_grid_list(my_file)
    vox_ind, vols = get_dend_indices(my_grid, region=region)
    specie_list = get_all_species(my_file, output)
    population = get_populations(my_file, trial, output)
    specie_idx = []
    for sp in specie:
        specie_idx.append(specie_list.index(sp))
    voxel_list = sorted(vox_ind.keys())
    how_many_voxels = len(voxel_list)
    out = np.zeros((population.shape[0], how_many_voxels))
    for i, key in enumerate(voxel_list):
        volume = vols[key]
        for idx in specie_idx:
            h = population[:, vox_ind[key], idx].sum(axis=1)
            out[:, i] += nano_molarity(h, volume)
    return out, voxel_list
        


if __name__ == '__main__':
    specie_list = ["Ca"]
    specie = "Ca"
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
        conc_dict = {}
        time_dict = {}
        for trial in my_file.keys():
            if trial == "model":
                continue
            conc, voxels = get_dynamics_in_region(my_file,
                                                  specie_list,
                                                  reg_list, trial, "__main__")
            conc_dict[trial] = conc
            time = get_times(my_file, trial, "__main__")
            time_dict[trial] = time
        vmax = 0
        vmin = 1000000000

        for key in conc_dict:
            new_max = conc_dict[key].max()
            new_min =  conc_dict[key].min()
            if new_max > vmax:
                vmax = new_max
            if new_min < vmin:
                vmin = new_min
        for i, key in enumerate(conc_dict):
            fig, ax = plt.subplots(1, 1)
            time = time_dict[key]
            im = ax.imshow(conc_dict[key].T, aspect="auto",
                           interpolation="none",
                           origin="lower", extent = [time[0]*1e-3,
                                                     time[-1]*1e-3,
                                                     voxels[0],
                                                     voxels[-1]],
                           cmap=plt.get_cmap("Reds"), vmin=vmin, vmax=vmax)
            fig.colorbar(im)
            ax.set_title("%s trial %d %s" % (fname, i, specie))
    plt.show()
                          
