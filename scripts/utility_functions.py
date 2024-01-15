#!/usr/bin/env python
import os
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import matplotlib.pyplot as plt

NA = Avogadro*1e-23
spine = ['PSD', 'head', 'neck']
t_init = 3000
window = 50

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


def nano_molarity(N, V):
    return 10 * N / V / NA


def pico_sd(N, S):
    return 10 * N / S / NA


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


def get_all_anchored_species(root):
    all_species = []
    for son in root:
        if son.tag.endswith('ReactionScheme'):
            for grandson in son:
                if grandson.tag.endswith('Specie'):
                    if not float(grandson.get("kdiff")):
                        all_species.append(grandson.get('id'))
    return list(set(all_species))


def get_all_submembrane_species(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    all_anchored_species = get_all_anchored_species(root)
    anchored = []
    for son in root:
        if son.tag.endswith('InitialConditions'):
            for grandson in son:
                if grandson.tag.endswith("SurfaceDensitySet"):
                    for grandgrandson in grandson:
                        name = grandgrandson.get("specieID")
                        if name in all_anchored_species:
                            anchored.append(name)
    return list(set(anchored))

def get_output_regions(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    outputs = {}
    for son in root:
        if son.tag.endswith('OutputScheme'):
            for grandson in son:
                outputs[grandson.get("filename")] = grandson.get("region")
    return outputs
        
def get_key(cell):
    if cell[18]:
        return cell[15].decode('utf-8') + '_' + cell[18].decode('utf-8')
    return cell[15].decode('utf-8')

def region_volumes(my_file):
    if isinstance(my_file, str):
        my_file = h5py.File(my_file)
    grid_list = get_grid_list(my_file)
    regions = get_regions(my_file)
    volumes = {}
    for region in regions:
        volumes[region] = 0
    for cell in grid_list:
        key = get_key(cell)
        volumes[key] += float(cell[12])

    return volumes


def sum_volume(my_file, region_list):
    grid_list = get_grid_list(my_file)
    vol_sum = 0
    volumes = region_volumes(my_file)
    for region in region_list:
        if region in volumes:
            vol_sum += volumes[region]
    return vol_sum


def sum_indices(my_file, region_list):
    reg_indices = get_region_indices(my_file)
    sum_indices = []
    for region in region_list:
        if region in reg_indices:
            sum_indices += reg_indices[region]
    return sum_indices


def region_surface(grid_list, direction=0):
    submembrane_regions = []
    submembrane_regions_dict = {}
    for i, cell in enumerate(grid_list):
        if cell[17] == b'submembrane':
            new_name = cell[15].decode('utf-8')
            if  new_name not in submembrane_regions:
                submembrane_regions.append(new_name)
                submembrane_regions_dict[new_name] = []
            submembrane_regions_dict[new_name].append(i)
    surface = {}
    for key in submembrane_regions_dict:
        surface[key] = 0
        for cell_idx in submembrane_regions_dict[key]:
            if direction == 0:
                depth = grid_list[cell_idx][13]
                width = abs(grid_list[cell_idx][0] - grid_list[cell_idx][3])
                surface[key] += depth * width
            else:
                print('Unimplemented direction', direction)

    return surface


def get_region_indices(my_file):
    grid_list = get_grid_list(my_file)
    region_ind = {}
    for idx, cell in enumerate(grid_list):
        key = get_key(cell)
        if key not in region_ind:
            region_ind[key] = []
        region_ind[key].append(idx)
    return region_ind

def get_spines(regions):
    out = {}
    for region in regions:
        try:
            end = region.split("_")[1]
        except IndexError:
            continue
        if end == "":
            continue

        if end in out:
            out[end].append(region)
        else:
            out[end] = [region]
    return out


def get_regions(my_file):
    grid_list = get_grid_list(my_file)
    return sorted(list(set([get_key(grid) for grid in grid_list])))


def get_concentrations_region_list(my_file, my_list, trial, out):

    grid_list = get_grid_list(my_file)
    species = get_all_species(my_file)
    idxs = sum_indices(my_file, my_list)
    vol = sum_volume(my_file, my_list)
    data = get_populations(my_file, trial=trial, output=out)
    numbers = data[:, idxs, :].sum(axis=1)
    return nano_molarity(numbers, vol)


def get_concentrations(my_file, trial, out):
    grid_list = get_grid_list(my_file)
    data = get_populations(my_file, trial=trial, output=out)
    species = get_all_species(my_file, output=out)
    regions = get_regions(my_file)
    submembrane_species = get_all_submembrane_species(my_file)
    volume_dict = region_volumes(my_file)
    surface_dict = region_surface(grid_list)
    concentrations = np.zeros((data.shape[0], len(regions), len(species)))
    numbers = np.zeros_like(concentrations)
    region_indices = get_region_indices(my_file)

    for i, reg in enumerate(regions):
        # get numbers
        numbers[:, i, :] = data[:, region_indices[reg], :].sum(axis=1)
        if reg in surface_dict:
            for j, specie in enumerate(species):
                if specie in submembrane_species:
                    concentrations[:, i, j] = pico_sd(numbers[:, i, j],
                                                      surface_dict[reg])
                else:
                    concentrations[:, i, j] = nano_molarity(numbers[:, i, j],
                                                            volume_dict[reg])
        else:
            concentrations[:, i, :] = nano_molarity(numbers[:, i, :],
                                                    volume_dict[reg])

        
    return concentrations


def save_single_file(times, concentrations, species, fname):
    header = 'time'
    for specie in species:
        header += ' ' + specie
    what_to_save = np.zeros((concentrations.shape[0], len(species) + 1))
    what_to_save[:, 0] = times[:concentrations.shape[0]]
    what_to_save[:, 1:] = concentrations
    print(fname)
    np.savetxt(fname, what_to_save, header=header, comments='')


def save_concentrations(my_file, fname_base, output, trial='trial0'):
    regions = get_regions(my_file)
    times = get_times(my_file, trial=trial, output=output)
    species = get_all_species(my_file, output=output)
    concentrations = get_concentrations(my_file, trial, output)
    if output == '__main__':
        add = ''
    else:
        add = output + '_'
    for i, region in enumerate(regions):
        fname = '%s_%s%s_%s.txt' % (fname_base, add, trial, region)
        save_single_file(times, concentrations[:, i, :], species, fname)
    if len(regions) > 1:
        totals = get_concentrations_region_list(my_file, regions, trial, output)
        save_single_file(times, totals, species,
                         '%s_%s%s_%s.txt' % (fname_base, add, trial, 'total'))
        spines_dict = get_spines(regions)
        for spine_name in spines_dict.keys():
            spine_reg = spines_dict[spine_name]
            spine = get_concentrations_region_list(my_file, spine_reg,
                                                   trial, output)
            save_single_file(times, spine, species,
                             '%s_%s%s_%s.txt' % (fname_base, add,
                                                 trial, spine_name))


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


def get_conc(fullname, specie_list, region_list, output_name):
    print(fullname)
    try:
        my_file = h5py.File(fullname)
    except FileNotFoundError:
        return
    conc_dict = {}
    time_dict = {}

    for trial in my_file.keys():
        if trial == "model":
            continue
        conc, voxels = get_dynamics_in_region(my_file,
                                                    specie_list,
                                                    region_list, trial,
                                                    output_name)
        conc_dict[trial] = conc
        time = get_times(my_file, trial, output_name)
        time_dict[trial] = time
    lmin = min([len(conc) for conc in conc_dict.values()])
    
    time_end = min([time[-1] for time in time_dict.values()])
    time_len = min([len(time) for time in time_dict.values()])
    time = np.linspace(0, time_end, time_len)
    shape2 = max([conc.shape[1] for conc in conc_dict.values()])
    conc_mean = np.zeros((lmin, shape2))
    for conc in conc_dict.values():
        conc_mean[:lmin, :] += conc[:lmin, :]
    conc_mean /= len(conc_dict)
    return voxels, time, conc_mean


def extract_max_delay(concentration, dt):
    mean = concentration[:,:int(t_init/dt)-1].mean()
    length = concentration.shape[0]
    distance = np.linspace(-length/2, length/2, length)
    branch = concentration[:, int(t_init/dt):].max(axis=1)
    delay = np.zeros_like(branch)
    for idx in range(51, 102, 1):
        try:
            if branch[idx] > 1.5*mean:
                delay[idx] = concentration[idx, int(t_init/dt):].argmax()*dt

                #   print(idx, delay[idx], branch[idx], 1.5*mean)
            else:
                break
        except ValueError:
            break
    for idx in range(50, -1, -1):
        try:
            if branch[idx] > 1.5*mean:
                delay[idx] = concentration[idx, int(t_init/dt):].argmax()*dt
                #   print(idx, delay[idx], branch[idx], 1.5*mean)
            else:
                break
        except ValueError:
            break
    
    return distance, branch, delay

def ca_wave_propagation_figs(directiories_list, descr, dend_dict,
                       what_species, region_list, output_name, color_list,
                       label_list, what_type, marker_list):

    fig1, ax1 = plt.subplots(2, len(dend_dict), figsize=(21, 10))
    for k, directory in enumerate(directiories_list):
        my_path = os.path.join("..", directory)
        if len(descr):
            description = descr[directory] 
        im_list = {}
        for i, key in enumerate(dend_dict.keys()):
            im_list[key] = []
            for j, fname in enumerate(dend_dict[key]):
                if len(descr):
                    new_fname = fname % description
                else:
                    new_fname = fname
                my_file = os.path.join(my_path, new_fname)
                try:
                    vox, times, conc_mean = get_conc(my_file, ["Ca"],
                                                     region_list, output_name)
                except TypeError:
                    print("No file", new_fname)
                    continue
                
                im_list[key].append(conc_mean.T)
                dt = times[1]-times[0]
            for j, conc in enumerate(im_list[key]):
                distance, branch, delay = extract_max_delay(conc, dt)
                if k % 2:
                    symbol = "o"
                else:
                    symbol = "d"

                ax1[0][i].plot(distance, branch, color_list[j],
                               marker=symbol,
                               label=label_list[key][j]+" "+
                               what_type[directory],
                               linestyle="", fillstyle=marker_list[directory])
                ax1[1][i].plot(distance, delay, color_list[j], marker=symbol,
                               label=label_list[key][j]+" "+
                               what_type[directory],
                               linestyle="", fillstyle=marker_list[directory])
                    
            ax1[0][0].set_ylabel("Calcium amplitude (nM)", fontsize=15)
            ax1[0][i].set_title("Injection %s" % key, fontsize=15)
            ax1[0][i].set_yscale("log")
            ax1[1][i].set_xlabel("Distance from stimulated site (um)",
                                 fontsize=15)
            ax1[1][0].set_ylabel("Ca wave delay (ms)", fontsize=15)
            
            
    
    ax1[0][2].legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    for i, axes in enumerate(ax1):
        ylim2 = max([max(ax.get_ylim()) for ax in axes])
        ylim1 = min([min(ax.get_ylim()) for ax in axes])
        for j, ax in enumerate(axes):
            ax.set_ylim([ylim1, ylim2])
            if j:
                ax.set_yticks([])
            if i<1:
                ax.set_xticks([])
            ax.tick_params(axis='both', which='major', labelsize=14)
    return fig1


def ca_wave_propagation_figs_bal_tubes(directiories_list,
                                       descr, dend_dict,
                                       what_species, region_list,
                                       output_name,
                                       color_list,
                                       label_list, what_type,
                                       marker_list):
    symbol = "d"
    fig1, ax1 = plt.subplots(2, len(dend_dict), figsize=(20, 10))
    for k, directory in enumerate(directiories_list):
        my_path = os.path.join("..", directory)
        if len(descr):
            description = descr[directory] 
        im_list = {}
        for i, key in enumerate(dend_dict.keys()):
            im_list[key] = []
            for j, fname in enumerate(dend_dict[key]):
                if len(descr):
                    new_fname = fname % description
                else:
                    new_fname = fname
                my_file = os.path.join(my_path, new_fname)
                try:
                    vox, times, conc_mean = get_conc(my_file, ["Ca"],
                                                     region_list, output_name)
                except TypeError:
                    print("No file", new_fname)
                    continue

                
                im_list[key].append(conc_mean.T)
                dt = times[1] - times[0]
           
            for j, conc in enumerate(im_list[key]):
                distance, branch, delay = extract_max_delay(conc, dt)
                if j > 2:
                    fillstyle = "none"
                else:
                    fillstyle = "full"
                ax1[0][i].plot(distance, branch, color_list[j],
                               marker=symbol,
                               label=label_list[key][j]+" "+
                               what_type[directory],
                               linestyle="", fillstyle=fillstyle)
                ax1[1][i].plot(distance, delay, color_list[j], marker=symbol,
                               label=label_list[key][j]+" "+
                               what_type[directory],
                               linestyle="", fillstyle=fillstyle)
                    
            ax1[0][0].set_ylabel("Calcium amplitude (nM)", fontsize=15)
            ax1[0][i].set_title("Injection %s" % key, fontsize=15)
            ax1[0][i].set_yscale("log")
            ax1[1][i].set_xlabel("Distance from stimulated site (um)",
                                 fontsize=15)
            ax1[1][0].set_ylabel("Ca wave delay (ms)", fontsize=15)
            
            
    
    ax1[0][2].legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    for i, axes in enumerate(ax1):
        ylim2 = max([max(ax.get_ylim()) for ax in axes])
        ylim1 = min([min(ax.get_ylim()) for ax in axes])
        for j, ax in enumerate(axes):
            ax.set_ylim([ylim1, ylim2])
            if j:
                ax.set_yticks([])
            if i<1:
                ax.set_xticks([])
            ax.tick_params(axis='both', which='major', labelsize=14)
    return fig1



def ca_wave_propagation_figs_different_paradigms(directiories_list,
                                                 descr, paradigm_dict,
                                                 what_species, region_list,
                                                 output_name, color_list,
                                                 label_list, what_type,
                                                 marker_list):
    symbol = "d"
    fig1, ax1 = plt.subplots(2, len(paradigm_dict), figsize=(20, 10))
    for k, directory in enumerate(directiories_list):
        my_path = os.path.join("..", directory)
        if len(descr):
            description = descr[directory] 
        im_list = {}
        for i, key in enumerate(paradigm_dict.keys()):
            im_list[key] = []
            for j, fname in enumerate(paradigm_dict[key]):
                
                if len(descr):
                    new_fname = fname % description
                else:
                    new_fname = fname
                print(i, j, k, new_fname)
                my_file = os.path.join(my_path, new_fname)
                try:
                    vox, times, conc_mean = get_conc(my_file, ["Ca"],
                                                     region_list, output_name)
                except TypeError:
                    print("No file", new_fname)
                    continue

                
                im_list[key].append(conc_mean.T)
                dt = times[1] - times[0]
           
            for j, conc in enumerate(im_list[key]):
                
                distance, branch, delay = extract_max_delay(conc, dt)
                if k%2:
                    fillstyle = "none"
                else:
                    fillstyle = "full"
                ax1[0][i].plot(distance, branch, color_list[j],
                               marker=symbol,
                               label=label_list[key][j]+" "+
                               what_type[directory],
                               linestyle="", fillstyle=fillstyle)
                ax1[1][i].plot(distance, delay, color_list[j], marker=symbol,
                               label=label_list[key][j]+" "+
                               what_type[directory],
                               linestyle="", fillstyle=fillstyle)
                    
            ax1[0][0].set_ylabel("Calcium amplitude (nM)", fontsize=15)
            ax1[0][i].set_title("Injection duration %s ms" % key, fontsize=15)
            ax1[0][i].set_yscale("log")
            ax1[1][i].set_xlabel("Distance from stimulated site (um)",
                                 fontsize=15)
            ax1[1][0].set_ylabel("Ca wave delay (ms)", fontsize=15)
            
            
    try:
        ax1[0][2].legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    except IndexError:
        ax1[0][1].legend(loc='upper left', bbox_to_anchor=(1, 0.5))
        
    for i, axes in enumerate(ax1):
        ylim2 = max([max(ax.get_ylim()) for ax in axes])
        ylim1 = min([min(ax.get_ylim()) for ax in axes])
        for j, ax in enumerate(axes):
            ax.set_ylim([ylim1, ylim2])
            if j:
                ax.set_yticks([])
            if i<1:
                ax.set_xticks([])
            ax.tick_params(axis='both', which='major', labelsize=14)
    return fig1
