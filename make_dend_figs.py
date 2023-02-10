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
    reg_list = ["dend", "dend1", "dend2", "dend3", "dend4", "dend5", "dend6",
                "dend7", "dend8", "dend9", "dend10"]
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
                          
