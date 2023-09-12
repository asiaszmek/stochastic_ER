import sys
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import utility_functions as utils

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',  'tab:purple', 'tab:brown',
          'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', "b", "orange", "g",
          "r", "p" ]
stim_dend = "dend26"
dend_f = {
    "1.2":
    [
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_1.2_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_1.2_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_1.2_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_baloon_diam_1.2_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_baloon_diam_1.2_um_50_um_4000_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_1.2_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_1.2_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_1.2_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_tubes_diam_1.2_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_tubes_diam_1.2_um_50_um_4000_nM.h5",        
    ],
       "2.4":
    [
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_2.4_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_2.4_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_2.4_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_baloon_diam_2.4_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_baloon_diam_2.4_um_50_um_4000_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_2.4_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_2.4_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_2.4_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_tubes_diam_2.4_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_tubes_diam_2.4_um_50_um_4000_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_2.4_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_2.4_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_2.4_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_2.4_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_2.4_um_50_um_4000_nM.h5",    
        
    ],
       "6.0":
    [
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_6.0_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_6.0_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_baloon_diam_6.0_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_baloon_diam_6.0_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_baloon_diam_6.0_um_50_um_4000_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_6.0_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_6.0_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_tubes_diam_6.0_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_tubes_diam_6.0_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_tubes_diam_6.0_um_50_um_4000_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_6.0_um_50_um_350_nM.h5",
        "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_6.0_um_50_um_700_nM.h5",        
        "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_6.0_um_50_um_1050_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_6.0_um_50_um_2000_nM.h5",        
        # "model_RyR_3s_injection_simple_SERCA_nc_tubes_diam_6.0_um_50_um_4000_nM.h5",   
    ],
}
labels = {
    "1.2":
    [
        "dendritic membrane 350 nM",
        "dendritic membrane 700 nM",        
        "dendritic membrane 1050 nM",        
        "dendritic membrane 2000 nM",        
        "dendritic membrane 4000 nM",
        "RyR2 uniform  350 nM",
        "RyR2 uniform  700 nM",        
        "RyR2 uniform  1050 nM",        
        "RyR2 uniform  2000 nM",        
        "RyR2 uniform  4000 nM",        
    ],
       "2.4":
    [
        "RyR2 dendritic membrane 350 nM",
        "RyR2 dendritic membrane 700 nM",        
        "RyR2 dendritic membrane 1050 nM",        
        # "RyR2 dendritic membrane 2000 nM",        
        # "RyR2 dendritic membrane  4000 nM",
        "RyR2 uniform 350 nM",
        "RyR2 uniform 700 nM",        
        "RyR2 uniform 1050 nM",        
        # "RyR2 uniform 2000 nM",        
        # "RyR2 uniform 4000 nM",
        "RyR2 overexpressiom 350 nM",
        "RyR2 overexpression 700 nM",        
        "RyR2 overexpression 1050 nM",        
        # "RyR2 overexpression 2000 nM",        
        # "RyR2 overexpression 4000 nM",     
    ],
       "6.0":
    [
        "RyR2 dendritic membrane 350 nM",
        "RyR2 dendritic membrane 700 nM",        
        "RyR2 dendritic membrane 1050 nM",        
        # "RyR2 dendritic membrane 2000 nM",        
        # "RyR2 dendritic membrane 4000 nM",
        "RyR2 uniform 350 nM",
        "RyR2 uniform 700 nM",        
        "RyR2 uniform 1050 nM",        
        # "RyR2 uniform 2000 nM",        
        # "RyR2 uniform 4000 nM",
        "RyR2 overexpression 350 nM",
        "RyR2 overexpression 700 nM",        
        "RyR2 overexpression 1050 nM",        
        # "RyR2 overexpression 2000 nM",        
        # "RyR2 overexpression 4000 nM",   
    ],
}
    
NA = Avogadro*1e-23
specie_dict = {
    "Ca": ["Ca"],
    "CaOut": ["CaOut"],
    "CaER": ["CaER"],
    "RyRO": ["RyRO1", "RyRO2"],
    "STIM_CaER": ["STIM_2CaER"],
    "Orai": ["OraiSTIM_4", "Orai2STIM_4", "Orai3STIM_4"]
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

 
    im_list = {}
    for i, key in enumerate(dend_f.keys()):
        im_list[key] = []
        for j, fname in enumerate(dend_f[key]):
            try:
                my_file = h5py.File(fname, 'r')
            except FileNotFoundError:
                print(fname, " not found")
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
            delay = [(conc[max_idx_seg_side1, 3000:].argmax()
                      +conc[max_idx_seg_side2, 3000:].argmax())/2*dt]
            max_pre = np.mean(conc[:, :3000].max(axis=1))
            for idx in range(1, 51):
                distance.append(idx/2)
                peak = (conc[max_idx_seg_side1-idx, 3000:].max()
                               +conc[max_idx_seg_side2+idx, 3000:].max())/2
                branch.append(peak)
               
                if peak > 1:
                    delay.append((conc[max_idx_seg_side1-idx, 3000:].argmax()
                                  +conc[max_idx_seg_side2+idx, 3000:].argmax())/2*dt)
                    
                else:
                     delay.append(0)
                    
            ax1[0][i].plot(distance, branch, colors[j], marker = "d",
                           label=labels[key][j])
            ax1[1][i].plot(distance, delay, colors[j], marker = "d",
                           label=labels[key][j])
 
        # ax1[0][i].plot(distance, np.log10(np.ones_like(distance)*1e-7),
        #              "k", label = "100 nM")
        #ax1[0][i].set_xlabel("Distance from stimulated site (um)")
        ax1[0][0].set_ylabel("% basal calcium", fontsize=15)
        ax1[0][i].set_title("%s um diameter" % key, fontsize=15)
        
        ax1[1][i].set_xlabel("Distance from stimulated site (um)", fontsize=15)
        ax1[1][0].set_ylabel("Ca wave delay (ms)", fontsize=15)
       
    ax1[0][i].legend()
    

    for axes in ax1:
        ylim2 = max([max(ax.get_ylim()) for ax in axes])
        ylim1 = min([min(ax.get_ylim()) for ax in axes])
        for ax in axes:
            ax.set_ylim([ylim1, ylim2])
    fig1.savefig("Ca_wave_vs_distance_3s_injection.png",
                 bbox_inches=None, pad_inches=0.1)
    
    plt.show()
                          
