#!/usr/bin/env python
import os
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



NA = Avogadro*1e-23
spine = ['PSD', 'head', 'neck']
t_init = 3000
window = 50


def get_array(conc_dict, specie):
    mini  = min([conc_dict[specie][key].shape[-1]
                 for key in conc_dict[specie].keys()])

    no = len(conc_dict[specie].keys())
    voxels = conc_dict[specie]["trial0"].shape[0]
    out = np.zeros((no, voxels, mini))
    for i, trial in enumerate(conc_dict[specie].keys()):
        out[i] = conc_dict[specie][trial][:,:mini]
    return out
    

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


def get_dend_indices(grid, region=["dend"]):
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


def get_conc(fullname, specie_list, region_list, output):
    if isinstance(specie_list, str):
        specie_list = [specie_list]
    conc_dict = {}
    time_dict = {}
    for specie in specie_list:
        conc_dict[specie] = {}
    try:
        my_file = h5py.File(fullname)
        print(fullname)
    except FileNotFoundError:
        print("File not found", fullname)
        return conc_dict, time_dict
    
    for trial in my_file.keys():
        if trial == "model":
            continue
        try:
            for specie in specie_list:
                pop, voxel = get_dynamics_in_region(my_file, specie, region_list,
                                                    trial, output)
                conc_dict[specie][trial] = pop.T
            time = get_times(my_file, trial, output)
            time_dict[trial] = time
        except IOError:
            print("Something wrong with", fullname)
            break
    return conc_dict, time_dict

        
def extract_max_delay(conc_dict, dt):
    length = conc_dict["trial0"].shape[0]
    branch = np.zeros((len(conc_dict), length))
    delay = np.zeros((len(conc_dict), length))
    distance = np.linspace(-length/2, length/2, length)

    for i, concentration in enumerate(conc_dict.values()):
        maxi = concentration[:, :int(t_init/dt)-1].max()
        mean = concentration[:, :int(t_init/dt)-1].mean()
        std = concentration[:, :int(t_init/dt)-1].std()
        try:
            branch[i] = concentration[:, int(t_init/dt):].max(axis=1)
        except ValueError:
            continue
        for idx in range(51, 102, 1):
            try:
                if branch[i, idx] > maxi and  branch[i, idx] > mean+3*std:
                    delay[i, idx] = concentration[idx, int(t_init/dt):].argmax()*dt
                   
                else:
                    break
            except ValueError:
                break
        for idx in range(50, -1, -1):
            try:
                if branch[i, idx] > maxi and  branch[i, idx] > mean+3*std:
                    delay[i, idx] = concentration[idx, int(t_init/dt):].argmax()*dt
                 
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
                    times_dict, conc_dict = get_conc(my_file, ["Ca"],
                                                     region_list, output_name)
                except TypeError:
                    continue
                
                im_list[key].append(conc_dict)
                dt = times_dict["trial0"][1]-times_dict["trial0"][0]
            for j, conc in enumerate(im_list[key]):
                try:
                    distance, branch, delay = extract_max_delay(conc_dict["Ca"], dt)
                except TypeError:
                    continue
                
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
                    conc_dict, times_dict = get_conc(my_file,
                                                     ["Ca"],
                                                     region_list,
                                                     output_name)
                except TypeError:
                    continue

                
                im_list[key].append(conc_mean.T)
                dt = times_dict["trial0"][1] - times_dict["trial0"][0]
           
            for j, conc in enumerate(im_list[key]):
                try:
                    distance, branch, delay = extract_max_delay(conc_dict["Ca"], dt)
                except TypeError:
                    continue
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
                my_file = os.path.join(my_path, new_fname)
                try:
                    conc_dict, times_dict = get_conc(my_file, ["Ca"],
                                                     region_list, output_name)
                except TypeError:
                    continue

                
                im_list[key].append(conc_mean.T)
                dt = times_dict["trial0"][1] - times_dict["trial0"][0]
           
            for j, conc in enumerate(im_list[key]):
                try:
                    distance, branch, delay = extract_max_delay(conc_dict["Ca"], dt)
                except TypeError:
                    continue
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



def make_distance_fig(fname, directories_list, descr, dend_diam, stims,
                      what_species, region_list, output_name,
                      colors, types):
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    paradigm_dict = {
        0: "",
        1: "_3s_injection"}
    for k, d in enumerate(directories_list):
        my_path = os.path.join("..", d)
        add = descr[d]
        im_list = {}
        for l, inh in enumerate([""]):
            for j, diam in enumerate(dend_diam):
                y = []
                y_err = []
                x = []
                for i, stim in enumerate(stims):
                    new_fname = fname % (add, inh, diam, stim)
                    my_file = os.path.join(my_path, new_fname)
                    

                    conc_dict, times_dict = get_conc(my_file,
                                                         what_species,
                                                         region_list,
                                                         output_name)
                    try:
                        dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                    except KeyError:
                        continue
                    distance, branch, delay = extract_max_delay(conc_dict["Ca"],
                                                                    dt)

                    
                    full_delay = np.zeros((len(delay),))
                    for i, delay_1d in enumerate(delay): 
                        full_delay[i] = len(np.where(delay_1d>0)[0])/4
                    y.append(full_delay.mean())
                    y_err.append(full_delay.std()/len(full_delay)**0.5)
                    b_diam = float(diam)
                    x.append(np.mean(branch[:,50]+branch[:,51])/2000)
                print(x, y, y_err)
                if k == 0:
                    ax1.errorbar(x, y, yerr=y_err, color=colors[diam], marker="d",
                                label=types[d]+" diam "+diam,
                                linestyle="")

                if k == 1:
                    ax1.errorbar(x, y, yerr=y_err, color=colors[diam], marker="o",
                                label=types[d]+" diam "+diam,
                                linestyle="")
                if k == 2:
                    ax1.errorbar(x, y, yerr=y_err, color=colors[diam], marker="d",
                                label=types[d]+" diam "+diam, linestyle="",
                                fillstyle="none")
                if k == 3:
                    ax1.errorbar(x, y, yerr=y_err, color=colors[diam], marker="o",
                                label=types[d]+" diam "+diam, linestyle="",
                                fillstyle="none")
    ax1.set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    ax1.set_ylabel("Spatial extent [um]", fontsize=20)
    
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.legend(loc='lower left', bbox_to_anchor=(1, 0.5))
    return fig1


def make_distance_fig_2_4(fname, directories, descr, dend_diam,
                          stims, what_species, organization,
                          dur_dict, reg_list, output_name, 
                          colors, types, marker):
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))


    for d in directories:
        my_path = os.path.join("..", d)
        add = descr[d]
        for k, org in  enumerate(organization):
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (add, inh, org, diam, stim)
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
               
                
                        try:
                            distance, branch, delay = extract_max_delay(conc_dict["Ca"],
                                                                        dt)
                        except TypeError:
                            continue
                        

                    
                        full_delay = np.zeros((len(delay),))
                        for i, delay_1d in enumerate(delay): 
                            full_delay[i] = len(np.where(delay_1d>0)[0])/4
                        y.append(full_delay.mean())
                        y_err.append(full_delay.std()/len(full_delay)**0.5)
                        b_diam = float(diam)
                        x.append(np.mean(branch[:,50]+branch[:,51])/2000)
                    print(x, y, y_err)
                    if k % 2:
                        ax1.errorbar(x, y, yerr=y_err,
                                     color=colors[diam], marker=marker[inh],
                                     label=types[org]+" diam "+ diam + dur_dict[inh],
                                     linestyle="", fillstyle="full")
                    else:
                        ax1.errorbar(x, y, yerr=y_err, color=colors[diam],
                                     marker=marker[inh],
                                     label=types[org]+" diam "+diam
                                     + dur_dict[inh],
                                     linestyle="", fillstyle="none")
    ax1.set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    ax1.set_ylabel("Spatial extent [um]", fontsize=20)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return fig1



def make_distance_fig_aging(directories, descr, dend_diam,
                            stims, what_species, organization,
                            dur_dict, reg_list, output_name, 
                            colors, types, marker):
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    for k, d in enumerate(directories):
        my_path = os.path.join("..", d)
        add = descr[d]
        fname = directories[d]
        for org in organization:
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (add, inh, org, diam, stim)
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
               
                
                        try:
                            distance, branch, delay = extract_max_delay(conc_dict["Ca"],
                                                                        dt)
                        except TypeError:
                            continue
                        full_delay = np.zeros((len(delay),))
                        for i, delay_1d in enumerate(delay): 
                            full_delay[i] = len(np.where(delay_1d>0)[0])/4
                        y.append(full_delay.mean())
                        y_err.append(full_delay.std()/len(full_delay)**0.5)
                        b_diam = float(diam)
                        x.append(np.mean(branch[:,50]+branch[:,51])/2000)
                    print(x, y, y_err)
                    if not len(y):
                        continue
                    if not k % 2:
                        ax1.errorbar(x, y, yerr=y_err, color=colors[diam], marker=marker[inh],
                                     label=types[d]+" diam "+diam + dur_dict[inh],
                                     linestyle="", fillstyle="full")
                    else:
                        ax1.errorbar(x, y, yerr=y_err, color=colors[diam], marker=marker[inh],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[inh],
                                     linestyle="", fillstyle="none")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    ax1.set_ylabel("Spatial extent [um]", fontsize=20)

    return fig1


def make_distance_fig_aging_CaER(directories,  dend_diam,
                                 stims, what_species, organization,
                                 dur_dict, reg_list, output_name, 
                                 colors, types):
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    fig2, ax2 = plt.subplots(1, 1, figsize=(5, 5))
    stim_labels = {
        "": " 40 ms",
        "_3s_injection": " 3 ms"
    }
    marker = {
        "": "d",
        "_3s_injection": "o"
                                 
        }
    for k, d in enumerate(directories):
        my_path = os.path.join("..", d)
        fname = directories[d]
        for stim_type in ["", "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    x_err = []
                    y_ER = []
                    y_ER_err = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (stim_type, inh, diam, stim)
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca", "CaER"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
               
                
                        try:
                            distance_ca, branch_ca, delay_ca = extract_max_delay(conc_dict["Ca"],
                                                                        dt)
                        except TypeError:
                            continue
                
                       
                        full_delay = np.zeros((len(delay_ca),))
                        for i, delay_1d in enumerate(delay_ca): 
                            full_delay[i] = len(np.where(delay_1d>0)[0])/4
                        x.append(full_delay.mean())
                        x_err.append(full_delay.std()/len(full_delay)**.5)
                        full_dip = np.zeros((len(delay_ca),))
                        for i, key in enumerate(conc_dict["CaER"].keys()): 
                            full_dip[i] = (min(conc_dict["CaER"][key][50, int(t_init/dt):])
                                           +min(conc_dict["CaER"][key][51, int(t_init/dt):]))/2000*4 #  uM
                          
                        y.append(full_dip.mean())
                        y_err.append(full_dip.std()/len(full_dip)**0.5)
                        b_diam = float(diam)
                        myfile = h5py.File(my_file)
                        my_grid = get_grid_list(myfile)
                        vox_ind, vols = get_dend_indices(my_grid,
                                                         region=reg_list)
                        volume = sum(vols)
                        y_ER.append(y[-1]*volume*4*6.022)
                        y_ER_err.append(y_err[-1]*volume*4*6.022)
                    print(x, x_err, y_ER, y_ER_err, y, y_err)
                    if not len(y):
                        continue
                    if not k % 2:
                        ax1.errorbar(y, x, xerr=y_err, yerr=x_err,
                                     color=colors[diam],
                                     marker=marker[stim_type],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[inh]
                                     +stim_labels[stim_type],
                                     linestyle="", fillstyle="full")
                        ax2.errorbar(y_ER, x, xerr=y_ER_err, yerr=x_err,
                                     color=colors[diam],
                                     marker=marker[stim_type],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[inh]
                                     +stim_labels[stim_type],
                                     linestyle="", fillstyle="full")

                        print(types[d]+" diam "+diam + dur_dict[inh], k, "full")
                    else:
                        ax1.errorbar(y, x, yerr=x_err, xerr=y_err,
                                     color=colors[diam],
                                     marker=marker[stim_type],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[inh]
                                     +stim_labels[stim_type],
                                     linestyle="", fillstyle="none")
                        ax2.errorbar(y_ER, x, yerr=x_err, xerr=y_ER_err,
                                     color=colors[diam],
                                     marker=marker[stim_type],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[inh]
                                     +stim_labels[stim_type],
                                     linestyle="", fillstyle="none")

                        print(types[d]+" diam "+diam + dur_dict[inh],k, "none")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_xlabel("min Ca in the ER [uM]", fontsize=20)
    ax1.set_ylabel("Spatial extent [um]", fontsize=20)
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_xlabel("min Ca molecules in the ER",
                   fontsize=20)
    ax2.set_ylabel("Spatial extent [um]", fontsize=20)
    

    return fig1, fig2



def fit_exp(time, ca_conc, dt):
    time = time[int(t_init/dt):] - t_init
    ca_conc = ca_conc[int(t_init/dt):] - ca_conc[:int(t_init/dt)].mean()
    max_idx = ca_conc.argmax()
    duration = 2000

    ca_conc_log = ca_conc[max_idx:duration+max_idx]
    new_time = time[max_idx:duration+max_idx]-time[max_idx]
                         
    popt, pcov = curve_fit(lambda t, a, b, c: a*np.exp(-t/b)+c,
                           new_time, ca_conc_log)
    return popt[1]

def make_decay_constant_fig(directories,  dend_diam,
                            stims, what_species, organization,
                            dur_dict, output_name, 
                            colors, types):
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))

    stim_labels = {
        "": " 40 ms",
        "_3s_injection": " 3 ms"
    }
    marker = {
        "": "d",
        "_3s_injection": "o"
                                 
        }
    reg_list = ["dend26"]
    for k, d in enumerate(directories):
        my_path = os.path.join("..", d)
        fname = directories[d]
        for stim_type in [""]:#, "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (stim_type, inh, diam, stim)
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
                        ca_means = np.zeros((len(conc_dict["Ca"].keys())))
                        t_decays1 = np.zeros((len(conc_dict["Ca"].keys())))
                        for i, trial in enumerate(conc_dict["Ca"].keys()):
                            time = times_dict[trial]
                            ca = conc_dict["Ca"][trial].mean(axis=0)
                            t1 = fit_exp(time, ca, dt)
                            t_decays1[i] = t1
                            ca_means[i] = ca.max()/1000
                        x.append(ca_means.mean())
                        y.append(t_decays1.mean())
                        y_err.append(t_decays1.std()/len(t_decays1)**0.5)
                        
                        
                    print(x, y, y_err)
                    if not len(y):
                        continue
                    if not k % 2:
                        ax1.errorbar(x, y,  yerr=y_err,
                                     color=colors[diam],
                                     marker=marker[stim_type],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[stim_type]
                                     +stim_labels[stim_type],
                                     linestyle="", fillstyle="full")
                    else:
                        ax1.errorbar(x, y, yerr=y_err,
                                     color=colors[diam],
                                     marker=marker[stim_type],
                                     label=types[d]+" diam "+diam
                                     + dur_dict[stim_type]
                                     +stim_labels[stim_type],
                                     linestyle="", fillstyle="none")

    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    ax1.set_ylabel("Temporal decay constant [m sec]", fontsize=20)
    return fig1



def make_decay_constant_fig_sep_dends(directories,  dend_diam,
                                      stims, what_species, organization,
                                      dur_dict, output_name, 
                                      colors, types):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(15, 5))
    if len(dend_diam) == 1:
        ax1 = [ax1]
    
    marker = {
        "": "d",
        "_3s_injection": "o"
                                 
        }
    reg_list = ["dend26"]
    for k, d in enumerate(directories):
        my_path = os.path.join("..", d)
        fname = directories[d]
        for stim_type in [""]:  #  , "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (stim_type, inh, diam, stim)
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
                        ca_means = np.zeros((len(conc_dict["Ca"].keys())))
                        t_decays1 = np.zeros((len(conc_dict["Ca"].keys())))
                        for i, trial in enumerate(conc_dict["Ca"].keys()):
                            time = times_dict[trial]
                            ca = conc_dict["Ca"][trial].mean(axis=0)
                            try:
                                t1 = fit_exp(time, ca, dt)
                            except ValueError:
                                continue
                            t_decays1[i] = t1
                            ca_means[i] = ca.max()/1000
                        x.append(ca_means.mean())
                        y.append(t_decays1.mean())
                        y_err.append(t_decays1.std()/len(t_decays1)**0.5)
                        
                        
                    print(x, y, y_err)
                    if not len(y):
                        continue
                    
                    ax1[j].errorbar(x, y,  yerr=y_err,
                                    color=colors[d],
                                    marker=marker[stim_type],
                                    label=types[d]+ dur_dict[stim_type],
                                     linestyle="", fillstyle="full")
              
    ax1[-1].legend(loc=1)

    ax1[0].set_ylabel("Temporal decay constant [m sec]", fontsize=20)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])
    ax1[0].set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title("dend diam "+diam+  " um", fontsize=20)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
        
            
    return fig1


def make_spatial_specificity_fig_sep_dends(directories,  dend_diam,
                                           stims, what_species,
                                           dur_dict, output_name, 
                                           colors, types):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(15, 5))
    if len(dend_diam) == 1:
        ax1 = [ax1]
    stim_labels = {
        "": " 40 ms",
        "_3s_injection": " 3 ms"
    }
    marker = {
        "": "d",
        "_3s_injection": "o"
                                 
        }
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04",
                "dend05", "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
 
    for k, d in enumerate(directories):
        my_path = os.path.join("..", d)
        fname = directories[d]
        for stim_type in [""]:#, "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (stim_type, inh, diam, stim)
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
                        
                        try:
                            distance, branch, delay = extract_max_delay(conc_dict["Ca"],
                                                                        dt)
                        except TypeError:
                            continue
                        full_delay = np.zeros((len(delay),))
                        for i, delay_1d in enumerate(delay): 
                            full_delay[i] = len(np.where(delay_1d>0)[0])/4
                        y.append(full_delay.mean())
                        y_err.append(full_delay.std()/len(full_delay)**0.5)
                        b_diam = float(diam)
                        x.append(np.mean(branch[:,50]+branch[:,51])/2000)
                    print(x, y, y_err)
                    if not len(y):
                        continue

                    ax1[j].errorbar(x, y,  yerr=y_err,
                                    color=colors[d],
                                    marker=marker[stim_type],
                                    label=types[d]+ dur_dict[stim_type]
                                    +stim_labels[stim_type],
                                    linestyle="", fillstyle="full")


    ax1[-1].legend(loc=1)
    #ax1[0].legend(loc=1)
    ax1[0].set_ylabel("Spatial extent [um]", fontsize=20)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])
    ax1[0].set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title("dend diam "+diam+  " um", fontsize=20)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
        
            
    return fig1


def make_decay_constant_fig_ctrl(fname, directory,  dend_diam,
                                 stims, organization,
                                 dur_dict, output_name, 
                                 colors, types):
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    if len(dend_diam) == 1:
        ax1 = [ax1]
    stim_labels = {
        "": " 40 ms",
        "_3s_injection": " 3 ms"
    }
    marker = {
        "": "d",
        "_3s_injection": "o"
                                 
        }
    base = "dend"
    reg_list = ["dend26"]

    d = directory
    for k, org in enumerate(organization):
        my_path = os.path.join("..", d)
        for stim_type in [""]:
            for j, diam in enumerate(dend_diam):
                
                y = []
                y_err = []
                x = []
                for i, stim in enumerate(stims):
                    new_fname = fname % (stim_type, org, diam, stim)
                    my_file = os.path.join(my_path, new_fname)
                    try:
                        conc_dict, times_dict = get_conc(my_file,
                                                         ["Ca"],
                                                         reg_list,
                                                         output_name)
                    except TypeError:
                        continue
                    try:
                        dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                    except KeyError:
                        continue
                    ca_means = np.zeros((len(conc_dict["Ca"].keys())))
                    t_decays1 = np.zeros((len(conc_dict["Ca"].keys())))
                    for i, trial in enumerate(conc_dict["Ca"].keys()):
                        time = times_dict[trial]
                        ca = conc_dict["Ca"][trial].mean(axis=0)
                        t1 = fit_exp(time, ca, dt)
                        t_decays1[i] = t1
                        ca_means[i] = ca.max()/1000
                    x.append(ca_means.mean())
                    y.append(t_decays1.mean())
                    y_err.append(t_decays1.std()/len(t_decays1)**0.5)
                        
                print(x, y, y_err)
                if not len(y):
                    continue
                if not k % 2:
                    ax1.errorbar(x, y,  yerr=y_err,
                                 color=colors[d][diam],
                                 marker=marker[stim_type],
                                 label=types[org]+ dur_dict[stim_type]
                                 +stim_labels[stim_type],
                                 linestyle="", fillstyle="full")
                else:
                    ax1.errorbar(x, y, yerr=y_err,
                                 color=colors[d][diam],
                                 marker=marker[stim_type],
                                 label=types[org]
                                 + dur_dict[stim_type]
                                 +stim_labels[stim_type],
                                 linestyle="", fillstyle="none")

    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Temporal decay constant [m sec]", fontsize=20)
    ax1.set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    return fig1



def make_spatiotemporal_specificity_fig_sep_dends(directories,  dend_diam,
                                                  stims, what_species,
                                                  dur_dict, output_name, 
                                                  colors, types):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(15, 5))
    if len(dend_diam) == 1:
        ax1 = [ax1]
    stim_labels = {
        "": " 40 ms",
        "_3s_injection": " 3 ms"
    }
    marker = {
        "": "d",
        "_3s_injection": "o"
                                 
        }
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04",
                "dend05", "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
 
    for k, d in enumerate(directories):
        my_path = os.path.join("..", d)
        fname = directories[d]
        
        for stim_type in ["", "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                for inh in what_species:
                    y = []
                    y_err = []
                    x = []
                    x_err = []
                    for i, stim in enumerate(stims):
                        new_fname = fname % (stim_type, inh, diam, stim)
                        if new_fname == "model_RyR_simple_SERCA_tubes_diam_1.2_um_50_um_0175_nM.h5":
                            continue
                        my_file = os.path.join(my_path, new_fname)
                        try:
                            conc_dict, times_dict = get_conc(my_file,
                                                             ["Ca"],
                                                             reg_list,
                                                             output_name)
                        except TypeError:
                            continue
                        try:
                            dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        except KeyError:
                            continue
                        
                        try:
                            distance, branch, delay = extract_max_delay(conc_dict["Ca"],
                                                                        dt)
                        except TypeError:
                            continue
                        full_delay = np.zeros((len(delay),))
                        for i, delay_1d in enumerate(delay): 
                            full_delay[i] = len(np.where(delay_1d>0)[0])/4
                        x.append(full_delay.mean())
                        x_err.append(full_delay.std()/len(full_delay)**0.5)
                        
                        ca_means = np.zeros((len(conc_dict["Ca"].keys())))
                        t_decays1 = np.zeros((len(conc_dict["Ca"].keys())))
                        for i, trial in enumerate(conc_dict["Ca"].keys()):
                            time = times_dict[trial]
                            ca = conc_dict["Ca"][trial][50:52,:].mean(axis=0)
                            try:
                                t1 = fit_exp(time, ca, dt)
                            except ValueError:
                                break
                            t_decays1[i] = t1
                            ca_means[i] = ca.max()/1000
                        y.append(t_decays1.mean())
                        y_err.append(t_decays1.std()/len(t_decays1)**.5)
                    print(x,x_err, y, y_err)
                    if not len(y):
                        continue
                    
                    ax1[j].errorbar(x, y, xerr=x_err,
                                    yerr=y_err,
                                    color=colors[d],
                                    marker=marker[stim_type],
                                    label=types[d]+ dur_dict[stim_type]
                                    +stim_labels[stim_type],
                                    linestyle="", fillstyle="full")
                    

    ax1[-1].legend(loc=1)
    
    ax1[0].set_xlabel("Spatial extent [um]", fontsize=20)
    miniy = min([min(x.get_ylim()) for x in ax1])
    maxiy = max([max(x.get_ylim()) for x in ax1])
    minix = min([min(x.get_xlim()) for x in ax1])
    maxix = max([max(x.get_xlim()) for x in ax1])
    
    ax1[0].set_ylabel("Temporal decay constant [m sec]", fontsize=20)
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title("dend diam "+diam+  " um", fontsize=20)
        ax1[i].set_ylim([miniy, maxiy])
        ax1[i].set_xlim([minix, maxix])
        if i:
            ax1[i].set_yticks([])
        
            
    return fig1
