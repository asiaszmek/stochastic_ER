import sys
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import sys
from scipy.constants import Avogadro

NA = Avogadro*1e-23
specie_dict = {
    "Ca": ["Ca"],
    "CaOut": ["CaOut"],
    "CaER": ["CaER"],
    "RyRO": ["RyRO1", "RyRO2"],
    "STIM_CaER": ["STIM_2CaER"],
    "Orai": ["OraiSTIM_4", "Orai2STIM_4", "Orai3STIM_4"],
    "Fura": ["Fura2Ca"]
}
multiplier = {
    "Ca": 1,
    "CaOut": 1,
    "CaER": 1,
    "RyRO1": 1,
    "RyRO2": 1,
    "STIM_2CaER": 1,
    "OraiSTIM_4": 1,
    "Orai2STIM_4": 2,
    "Orai3STIM_4": 3,
    "Fura2Ca": 1,
}
def Parser():
    parser = argparse.ArgumentParser(description='Generate figs of avg conc')
    parser.add_argument('input', nargs='+',
                        help='input h5 files')
    parser.add_argument('--species', default="Ca",
                        help='Ca, RyRO, CaER, CaOut, RyRO, Fura')

    return parser

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
            h = population[:, vox_ind[key], idx].sum(axis=1)*multiplier[specie_list[idx]]
            out[:, i] += nano_molarity(h, volume)
    return out, voxel_list
        


if __name__ == '__main__':
    fnames = []
    args = Parser().parse_args()
    for name in args.input:
        if name.endswith("h5"):
            fnames.append(name)
    if not fnames:
        sys.exit('Do specify at least one totals filename')
    chosen_specie = args.species
    if chosen_specie in ["Ca", "CaER", "CaOut", "RyRO", "Fura"]:
        output_name = "all"
    elif  chosen_specie in ["STIM_CaER", "Orai"]:
        output_name = "RyR_Orai"

    try:
        specie_list = specie_dict[chosen_specie]
    except AttributeError:
        sys.exit("Unnkown specie %s" % s)
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
    fig1, ax1 = plt.subplots(1, len(fnames))
    im_list = []
    if len(fnames) == 1:
        ax1 = [ax1]
    for j, fname in enumerate(fnames):
        my_file = h5py.File(fname, 'r')
        conc_dict = {}
        time_dict = {}
        
        for trial in my_file.keys():
            if trial == "model":
                continue
            conc, voxels = get_dynamics_in_region(my_file,
                                                  specie_list,
                                                  reg_list, trial, output_name)
            conc_dict[trial] = conc
            time = get_times(my_file, trial, output_name)
            time_dict[trial] = time
        vmax = max([conc.max() for conc in conc_dict.values()])
        vmin = min([conc.min() for conc in conc_dict.values()])
       
        lmin = min([len(conc) for conc in conc_dict.values()])
        
        shape2 = max([conc.shape[1] for conc in conc_dict.values()])
        conc_mean = np.zeros((lmin, shape2))
        for conc in conc_dict.values():
            conc_mean[:lmin, :] += conc[:lmin, :]
        conc_mean /= len(conc_dict)

        # for i, key in enumerate(conc_dict):
        #     fig, ax = plt.subplots(1, 1)
        #     time = time_dict[key]
        #     im = ax.imshow(conc_dict[key].T, aspect="auto",
        #                    interpolation="none",
        #                    origin="lower", extent = [time[0]*1e-3,
        #                                              time[-1]*1e-3,
        #                                              voxels[0],
        #                                              voxels[-1]],
        #                    cmap=plt.get_cmap("Reds"), vmin=vmin, vmax=vmax)
        #     fig.colorbar(im)
        #     ax.set_title("%s trial %d %s" % (fname, i, specie))
        
        if chosen_specie == "Ca":
            im_list.append(np.log10(1e-9*conc_mean.T))
        else:
            im_list.append(conc_mean.T)
    vmax = max([im.max() for im in im_list])
    vmin = min([im.min() for im in im_list])
    for j, x in enumerate(ax1):
        im = ax1[j].imshow(im_list[j], aspect="auto",
                          interpolation="none",
                          origin="lower", extent = [time[0]*1e-3,
                                                    time[-1]*1e-3,
                                                    voxels[0],
                                                    voxels[-1]],
                          cmap=plt.get_cmap("Reds"), vmin=vmin,
                          vmax=vmax)

        if "baloon" in fnames[j]:
            ax1[j].set_title("baloon ER mean %s" %  chosen_specie)
        elif "tubes" in fnames[j]:
            ax1[j].set_title("stacked ER mean %s" %  chosen_specie)
            
    fig1.colorbar(im)
    plt.show()
                          
