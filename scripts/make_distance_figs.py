import os
import sys
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import utility_functions as utils
from scipy.constants import Avogadro


t_init = 3000

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',  'tab:purple', 'tab:brown',
          'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', "b", "orange", "g",
          "r", "p" ]
stim_dend = "dend26"
directories = [
    "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
    "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE",
]

symbol = {
    "tubes": "d",
    "baloon": "o",
}


dend_f = {
    "350 nM":
    [
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_1.2_um_50_um_0350_nM.h5",
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_0350_nM.h5",
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_6.0_um_50_um_0350_nM.h5",
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_1.2_um_50_um_0350_nM.h5",
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_2.4_um_50_um_0350_nM.h5",
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_6.0_um_50_um_0350_nM.h5",


           
    ],
       "700 nM":
    [
      
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_1.2_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_6.0_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_1.2_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_2.4_um_50_um_0700_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_6.0_um_50_um_0700_nM.h5",        
    
       
    ],
       "1050 nM":
    [
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_1.2_um_50_um_1050_nM.h5",               
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_2.4_um_50_um_1050_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_6.0_um_50_um_1050_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_1.2_um_50_um_1050_nM.h5",               
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_2.4_um_50_um_1050_nM.h5",        
        "model_RyR2CaM_simple_SERCA_SOCE_baloon_diam_6.0_um_50_um_1050_nM.h5",        
      
         
    ],
}
labels = {
    "350 nM":
    [
        "RyR2 uniform  1.2 um",
        "RyR2 uniform  2.4 um",        
        "RyR2 uniform  6.0 um",
        "RyR2 dendritic membrane  1.2 um",
        "RyR2 dendritic membrane  2.4 um",        
        "RyR2 dendritic membrane  6.0 um",        

    ],
       "700 nM":
    [

        "RyR2 uniform  1.2 um",
        "RyR2 uniform  2.4 um",        
        "RyR2 uniform  6.0 um",
        "RyR2 dendritic membrane  1.2 um",
        "RyR2 dendritic membrane  2.4 um",        
        "RyR2 dendritic membrane  6.0 um",        


    ],
       "1050 nM":
    [
        "RyR2 uniform  1.2 um",
        "RyR2 uniform  2.4 um",        
        "RyR2 uniform  6.0 um",        
        "RyR2 dendritic membrane  1.2 um",
        "RyR2 dendritic membrane  2.4 um",        
        "RyR2 dendritic membrane  6.0 um",        

    ],
}
    
NA = Avogadro*1e-23
specie_dict = {
    "Ca": ["Ca"],
    "CaOut": ["CaOut"],
    "CaER": ["CaER"],
    "RyRO": ["RyRO1", "RyRO2", "RyRCaMO1", "RyRCaMO2"],
    "STIM_CaER": ["STIM_2CaER"],
    "Orai": ["OraiSTIM_4", "Orai2STIM_4", "Orai3STIM_4"]
}
multiplier = {
    "Ca": 1,
    "CaOut": 1,
    "CaER": 1,
    "RyRO1": 1,
    "RyRO2": 1,
    "RyRCaMO1": 1,
    "RyRCaMO2": 1,
    "STIM_2CaER": 1,
    "OraiSTIM_4": 1,
    "Orai2STIM_4": 2,
    "Orai3STIM_4": 3,
}

def Parser():
    parser = argparse.ArgumentParser(description='Generate figs of avg conc')
    parser.add_argument('--species', default="Ca",
                        help='Ca, RyRO, CaER, CaOut, RyR')

    return parser

      


if __name__ == '__main__':
    fnames = []
    args = Parser().parse_args()
    chosen_specie = args.species
    if chosen_specie in ["Ca", "CaER", "CaOut", "RyRO"]:
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
        
    fig1, ax1 = plt.subplots(2, len(dend_f), figsize=(20, 10))
    for k, directory in enumerate(directories):
        my_path = os.path.join("..", directory)
        im_list = {}
        for i, key in enumerate(dend_f.keys()):
            im_list[key] = []
            for j, fname in enumerate(dend_f[key]):
                if k:
                    new_fname = fname.split("SOCE_")[0]+fname.split("SOCE_")[-1]
                    my_file = os.path.join(my_path, new_fname)
                else:
                    my_file = os.path.join(my_path, fname)
                try:
                    my_file = h5py.File(my_file, 'r')
                except FileNotFoundError:
                    print(my_file, " not found")
                    continue
                conc_dict = {}
                time_dict = {}
                for trial in my_file.keys():
                    if trial == "model":
                        continue
                    conc, voxels = utils.get_dynamics_in_region(my_file,
                                                                specie_list,
                                                                reg_list, trial, output_name)
                    conc_dict[trial] = conc
                    time = utils.get_times(my_file, trial, output_name)
                    time_dict[trial] = time
                    dt = time[1]-time[0]


                lmin = min([len(conc) for conc in conc_dict.values()])
        
                shape2 = max([conc.shape[1] for conc in conc_dict.values()])
                conc_mean = np.zeros((lmin, shape2))
                for conc in conc_dict.values():
                    conc_mean[:lmin, :] += conc[:lmin, :]
                conc_mean /= len(conc_dict)
                conc_mean = (conc_mean - conc_mean[:2000].mean(axis=0))/conc_mean[:2000].mean(axis=0)
                im_list[key].append(conc_mean.T)#np.log10(1e-9*conc_mean.T))
           
            for j, conc in enumerate(im_list[key]):
            

                length = conc.shape[0]
                distance = [0]
                max_idx_seg_side1 = 50
                max_idx_seg_side2 = 51
            
                branch = [(conc[max_idx_seg_side1].max()
                           +conc[max_idx_seg_side2].max())/2]
                delay = [(conc[max_idx_seg_side1, int(t_init/dt):].argmax()
                          +conc[max_idx_seg_side2, int(t_init/dt):].argmax())/2*dt]
                max_pre = np.mean(conc[:, :int(t_init/dt)].max(axis=1))
                for idx in range(1, 51):
                    distance.append(idx/2)
                    peak = (conc[max_idx_seg_side1-idx, int(t_init/dt):].max()
                            +conc[max_idx_seg_side2+idx, int(t_init/dt):].max())/2
                    branch.append(peak)
               
                    if peak > 1:
                        delay.append((conc[max_idx_seg_side1-idx, int(t_init/dt):].argmax()
                                      +conc[max_idx_seg_side2+idx, int(t_init/dt):].argmax())/2*dt)
                    
                    else:
                        delay.append(0)
                if j > 2:
                    symbol = "o"
                else:
                    symbol = "d"
                if not k:
                    ax1[0][i].plot(distance, branch, colors[j], marker=symbol,
                                   label=labels[key][j]+" SOCE", linestyle="")
                    ax1[1][i].plot(distance, delay, colors[j], marker=symbol,
                                   label=labels[key][j]+" SOCE", linestyle="")
                if k:
                    ax1[0][i].plot(distance, branch, colors[j], marker=symbol,
                                   label=labels[key][j], linestyle="", fillstyle="none")
                    ax1[1][i].plot(distance, delay, colors[j], marker=symbol,
                                   label=labels[key][j], linestyle="", fillstyle="none")
                    
            ax1[0][0].set_ylabel("% basal calcium", fontsize=15)
            ax1[0][i].set_title("Injection %s" % key, fontsize=15)
            
            ax1[1][i].set_xlabel("Distance from stimulated site (um)", fontsize=15)
            ax1[1][0].set_ylabel("Ca wave delay (ms)", fontsize=15)
            
            ax1[0][i].legend()
    

    for axes in ax1:
        ylim2 = max([max(ax.get_ylim()) for ax in axes])
        ylim1 = min([min(ax.get_ylim()) for ax in axes])
        for ax in axes:
            ax.set_ylim([ylim1, ylim2])
    fig1.savefig("Ca_wave_vs_distance.png",
                 bbox_inches=None, pad_inches=0.1)
    
    plt.show()
                          
