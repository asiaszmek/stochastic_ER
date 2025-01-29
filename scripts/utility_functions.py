#!/usr/bin/env python
import os
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


limit = 2.5
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


def get_length(My_file):
    if isinstance(My_file, str):
        try:
            my_file = h5py.File(My_file)
        except FileNotFoundError:
            return
    else:
        my_file = My_file    
    grid = get_grid_list(my_file)
    return grid[-1][3] - grid[0][0]
    
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


def fit_distance(conc_dict, dt, t_init=3000, method="regular", length=51):
    decays = np.zeros((len(conc_dict), 1))
    shape = conc_dict["trial0"].shape[0]
    distance = np.linspace(-length/2, length/2, shape)
    branch = np.zeros((len(conc_dict), 1))
    duration = conc_dict["trial0"].shape[1]
    if shape % 2:
        start_1 = start_2 = shape//2 + 1
    else:
        start_1 = shape//2-1
        start_2 = shape//2

    for i, concentration in enumerate(conc_dict.values()):
        ca_conc = np.zeros((shape,))
        ca_conc_mean = concentration[:, :int(t_init/dt)].mean()
        new_beg = int(t_init/dt)
        indices = []
        for j in range(start_2, shape):
            try:
                new_idx = concentration[j, new_beg:].argmax()
            except ValueError:
                continue
            ca_conc[j] = concentration[j, new_beg+new_idx]
            if method == "regular":
                if ca_conc[j] > limit*ca_conc_mean:
                    if not len(indices):
                        indices.append(j)
                    elif j+1 in indices or j-1 in indices:
                        indices.append(j)
        new_beg = int(t_init/dt)         
        for j in range(start_1, -1, -1):
            
            try:
                new_idx = concentration[j, new_beg:].argmax()
            except ValueError:
                continue
            ca_conc[j] = concentration[j, new_beg+new_idx]
            if method == "regular":
                if ca_conc[j] > limit*ca_conc_mean:
                    if not len(indices):
                        indices.append(j)
                    elif j+1 in indices or j-1 in indices:
                        indices.append(j)
        if method == "regular":
            dx = distance[1]-distance[0]
            decays[i] = len(indices)*dx/2
        else:
            try:
                popt, pcov = curve_fit(lambda t, a, b, c:
                                       a*np.exp(-abs(t)/b)+c,
                                       distance,
                                       ca_conc-(ca_conc[0]+ca_conc[-1])/2)
            except RuntimeError:
                continue
            if popt[1] < 0 or popt[1]> 10:
                continue
            decays[i] = popt[1]
        try:
            branch[i] = (concentration[start_1, int(t_init/dt):].max()
                         +concentration[start_2, int(t_init/dt):].max())/2
        except ValueError:
            continue
    return distance, branch, decays


def make_distance_fig_ratio(dir_list, directories_dict, dend_diam, stims,
                            what_species, region_list, output_name,
                            colors, types, method="regular"):
    # 0 -- denominator,
    # 1 -- numerator,
    # 2 -- denominator,
    # 3 -- numerator
    
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))


    paradigm_dict = {
        0: ""}
    l = 0
    for k, d in enumerate(dir_list):
        fname = directories_dict[d]
        my_path = os.path.join("..", d)
        
        im_list = {}
        if not k%2:
            res_den = {}
            res_den_error = {}
            x_den = {}
        if k%2:
            res_num = {}
            res_num_error = {}
            x_num = {}

        for j, diam in enumerate(dend_diam):
            if k%2 and diam not in res_num:
                res_num[diam] = []
                x_num[diam] = []
                res_num_error[diam]=[]
            elif not k%2 and diam not in res_den:
                res_den[diam] = []
                res_den_error[diam] = []
                x_den[diam] = []
  
            for i, stim in enumerate(stims):
                new_fname = fname % (diam, stim)
                my_file = os.path.join(my_path, new_fname)
                    

                conc_dict, times_dict = get_conc(my_file,
                                                 what_species,
                                                 region_list,
                                                 output_name)
                try:
                    dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                except KeyError:
                    continue
                length = get_length(my_file)
                try:
                    dist, branch, delay = fit_distance(conc_dict["Ca"],
                                                       dt, method=method,
                                                       length=length)
                except KeyError:
                    if k % 2:   #  numerator
                        res_num[diam].append([])
                        x_num[diam].append([])
                        res_num_error[diam].append([])
                if k % 2:   #  numerator
                    res_num[diam].append(delay.mean())
                    res_num_error[diam].append(delay.std()/len(delay)**0.5)
                    x_num[diam].append(branch)
                  
                else:
                    res_den[diam].append(delay.mean())
                    res_den_error[diam].append(delay.std()/len(delay)**0.5)
                    x_den[diam].append(branch)
            b_diam = float(diam)
            if k%2:
                out_num = np.array(res_num[diam])
                out_den = np.array(res_den[diam])
                err_num = np.array(res_num_error[diam])
                err_den = np.array(res_den_error[diam])

                ratio = out_num/out_den
                x_num[diam] = np.array(x_num[diam])
                x_den[diam] = np.array(x_den[diam])
                print(ratio)
                y = ratio
                x = (x_num[diam]+x_den[diam][:len(out_num)]).mean(axis=1)/2
                y_err = ((err_num/out_den)**2+(ratio/out_den*err_den)**2)**0.5
                x = x.reshape((len(x),))
                y = y.reshape((len(y),))
                y_err = y_err.reshape((len(y_err),))
            
            if k == 1:
                print(x, y, y_err)
                ax1[j].errorbar(x, y, yerr=y_err,  color=colors[diam],
                                marker="d",
                                label=types[d],
                                linestyle="", fillstyle="none")
            if k == 3:
                print(x, y, y_err)
                ax1[j].errorbar(x, y, yerr=y_err, color=colors[diam],
                                marker="o",
                                label=types[d], linestyle="",
                                fillstyle="none")
            if k == 5:
                print(x, y, y_err)
                ax1[j].errorbar(x, y, yerr=y_err, color=colors[diam],
                                marker="d",
                                label=types[d], linestyle="",
                                fillstyle="full")

            ax1[j].legend()
    ax1[0].set_xlabel("Peak Ca at stimulated site [uM]", fontsize=20)
    ax1[0].set_ylabel("Spatial extent ratio [um]", fontsize=20)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])

    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title("dend diam "+diam+  " um", fontsize=20)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
    return fig1


def make_injection_vs_min_CaER(directories,  dend_diam,
                               stims,  output_name, 
                               colors, types, marker, fillstyle):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04",
                "dend05", "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))

    for k, (d, fname) in enumerate(directories):
        my_path = os.path.join("..", d)
        for stim_type in [""]:    #, "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                x = []
                x_err = []
                y = []
                y_err = []
                for i, stim in enumerate(stims):
                    new_fname = fname % (stim_type, diam, stim)
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
                    branch = np.zeros((len(conc_dict["Ca"]), 1))
                    length = get_length(my_file)
                    dend_length = get_length(my_file)
                    shape = conc_dict["Ca"]["trial0"].shape[0]
                    if shape % 2:
                        start_1, start_2 = shape//2+1
                    else:
                        start_1 = shape//2-1
                        start_2 = shape//2
                    print(start_1, start_2)
                    min_mol = np.zeros((len(conc_dict["Ca"]), 1))
                    for i, key in enumerate(conc_dict["CaER"].keys()):
                        conc = conc_dict["CaER"][key][:, int(t_init/dt):].sum(axis=0)
                        conc_ca = conc_dict["Ca"][key]
                        branch[i] = (conc_dict["Ca"][key][start_1, int(t_init/dt):].max()+conc_dict["Ca"][key][start_2, int(t_init/dt):].max())/2000
                        min_mol[i] = min(conc)

                    b_diam = float(diam)
                    myfile = h5py.File(my_file)
                    my_grid = get_grid_list(myfile)
                    vox_ind, vols = get_dend_indices(my_grid,
                                                     region=reg_list)
                    volume = sum(vols)
                    x.append(np.mean(branch))
                    x_err.append(np.std(branch)/len(branch)**0.5)
                    y.append(min_mol.mean()*volume*4*6.022*1e-9)
                    y_err.append(min_mol.std()/len(min_mol)**.5*volume*4*6.022*1e-9)
                print(x,x_err, y, y_err)
                if not len(y):
                    continue
                ax1[j].errorbar(x, y, xerr=x_err, yerr=y_err,
                                color=colors[diam],
                                marker=marker[k],
                                label=types[k],
                                linestyle="", fillstyle=fillstyle[k])
                    
                ax1[j].legend(loc='lower left')

    ax1[0].set_ylabel("min Ca molecules in the ER [1e9]",
                   fontsize=20)
    ax1[0].set_xlabel("Peak Ca at stimulated site [um]", fontsize=20)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])
 
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title("dend diam "+diam+  " um", fontsize=20)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
    return fig1


def fit_exp(time, ca_conc, dt, duration=2000, t_init=3000, spatial=False):
    if not spatial:
        ca_conc_mean = ca_conc[:int(t_init/dt)].mean()
        ca_conc = ca_conc[int(t_init/dt):] - ca_conc_mean
        
        time = time[int(t_init/dt):] - t_init
        max_idx = ca_conc.argmax()
        ca_conc_log = ca_conc[max_idx:duration+max_idx]
        new_time = time[max_idx:duration+max_idx]-time[max_idx]
        try:
            popt, pcov = curve_fit(lambda t, a, b, c: a*np.exp(-t/b)+c,
                                   new_time, ca_conc_log)
        except RuntimeError:
            return 0
        return popt[1]
    popt, pcov = curve_fit(lambda t, a, b, c: a*np.exp(-abs(t)/b)+c,
                           time, ca_conc-(ca_conc[0]+ca_conc[-1])/2)

    return popt[1]
 

def make_decay_constant_fig_sep_dends(directories,  dend_diam,
                                      stims, output_name, 
                                      colors, types, marker, fillstyle):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(5*len(dend_diam), 5))
    if len(dend_diam) == 1:
        ax1 = [ax1]
    
    reg_list = ["dend26"]
    for k, (d, fname) in enumerate(directories):
        my_path = os.path.join("..", d)
        for stim_type in [""]:  #  , "_3s_injection"]:
            for j, diam in enumerate(dend_diam):
                y = []
                y_err = []
                x = []
                for i, stim in enumerate(stims):
                    new_fname = fname % (stim_type, diam, stim)
                    my_file = os.path.join(my_path, new_fname)
                    if diam == "1.2" and new_fname.startswith("model_RyR_simple_SERCA_SOCE_") and stim == "0175":
                        continue
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
                        if t1 > 0:
                            t_decays1[i] = t1
                            ca_means[i] = ca.max()/1000
                    x.append(ca_means.mean())
                    y.append(t_decays1.mean())
                    y_err.append(t_decays1.std()/len(t_decays1)**0.5)
                        
                        
                print(x, y, y_err)
                if not len(y):
                    continue

                ax1[j].errorbar(x, y,  yerr=y_err,
                                color=colors[diam],
                                marker=marker[k],
                                label=types[k],
                                fillstyle=fillstyle[k],
                                linestyle="")
                ax1[j].legend(loc="upper left")

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


def make_distance_fig_sep_dends(directories,  dend_diam,
                                stims, 
                                output_name, 
                                colors, types, marker, fillstyle,
                                method="regular"):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))
    if len(dend_diam) == 1:
        ax1 = [ax1]

    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04",
                "dend05", "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
 
    for k, (d, fname) in enumerate(directories):
        my_path = os.path.join("..", d)
        for stim_type in [""]:
            for j, diam in enumerate(dend_diam):
                y = []
                y_err = []
                x = []
                for i, stim in enumerate(stims):
                    new_fname = fname % (stim_type, diam, stim)
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
                    length = get_length(my_file)
                    try:
                        dist, branch, delay = fit_distance(conc_dict["Ca"],
                                                           dt,
                                                           method=method,
                                                           length=length)
                                                                    
                    except TypeError:
                        continue
                    y.append(delay.mean())
                    y_err.append(delay.std()/len(delay)**0.5)
                    b_diam = float(diam)
                    x.append(np.mean(branch)/1000)
                print(x, y, y_err)
                if not len(y):
                    continue

                ax1[j].errorbar(x, y,  yerr=y_err,
                                color=colors[diam],
                                fillstyle=fillstyle[k],
                                label=types[k], marker=marker[k],
                                linestyle="")
                if len(directories) == 4 and j == 2:
                    ax1[j].legend(loc="center right")
                else:
                    ax1[j].legend(loc="lower right")
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


# fig.suptitle(path)
# ax[0].imshow(new_conc,  aspect="auto",
#              interpolation="none",
#              origin="lower", extent = [0, 1000, 0.5, 51],
#              cmap=plt.get_cmap("Reds"))
# ax[1].imshow(edges,  aspect="auto",
#              interpolation="none",
#              origin="lower", extent = [0, 1000, 0.5, 51],
#              cmap=plt.get_cmap("Reds"))
# ax[2].imshow(new_conc,  aspect="auto",
#              interpolation="none",
#              origin="lower", extent = [0, 1000, 0.5, 51],
#              cmap=plt.get_cmap("Reds"))

# for clust in new_clusters:
#     #print(clust)
#     p0, p1 = clust
#     ax[2].plot((p0[1]*0.2, p1[1]*0.2), (p0[0]/2+1/2,
#                                         p0[0]/2+1/2),
#                linewidth=1)
#     ax[2].plot((p0[1]*0.2, p1[1]*0.2), (p1[0]/2+1/2,
#                                         p1[0]/2+1/2),
#                linewidth=1)
#     ax[2].plot((p0[1]*0.2, p0[1]*0.2), (p0[0]/2+1/2,
#                                         p1[0]/2+1/2),
#                linewidth=1)
#     ax[2].plot((p1[1]*0.2, p1[1]*0.2), (p0[0]/2+1/2,
#                                         p1[0]/2+1/2),
#                linewidth=1)

# plt.show()
