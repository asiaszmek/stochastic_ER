import os

import numpy as np
import pandas as pd
import utility_functions as utils



                  
colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}
directories = [
    [
        "Ca_wave_no_RyR_simple_SERCA_SOCE",
        "model_noRyR%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
        "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE",
        "model_RyR2CaM%s_simple_SERCA_SOCE_baloon_diam_%s_um_50_um_%s_nM.h5"
    ],
   
    [
        "Ca_wave_simple_SERCA_SOCE",
        "model_RyR%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_simple_SERCA_SOCE",
        "model_RyR%s_simple_SERCA_SOCE_baloon_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_normal_SERCA_aging",
        "model_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_normal_SERCA_aging",
        "model_aging%s_simple_SERCA_baloon_diam_%s_um_50_um_%s_nM.h5"
    ]
]
base = "dend"
reg_list = [base, "dend01", "dend02", "dend03", "dend04",
            "dend05", "dend06", "dend07", "dend08", "dend09",]
for i in range(10, 102, 1):
    reg_list.append("%s%d" %(base, i))
stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]


if __name__ == '__main__':
    types = [ "blocked RyR2", "uniform RyR2CaM", "RyR2CaM in EPJ","RyR2 no CaM uniform",
              "RyR2 no CaM in EPJ", "aging RyR2 uniform", "aging RyR2 in EPJ"]
    markers = ["o", "^", "o"]
    fillstyle = ["full", "full", "none"]
    stim_type = ""
    out_data = np.ndarray((len(directories), len(dend_diam), len(stims)))
    out_data_error = np.ndarray((len(directories), len(dend_diam), len(stims)))
    save_fname = "Spatial_specificity_%s_%s.csv"
    for j, diam in enumerate(dend_diam):
        for k, (d, fname) in enumerate(directories):
            my_path = os.path.join("..", d)
            for i, stim in enumerate(stims):
                new_fname = fname % (stim_type, diam, stim)
                my_file = os.path.join(my_path, new_fname)
                try:
                    conc_dict, times_dict = utils.get_conc(my_file,
                                                           ["Ca"],
                                                           reg_list,
                                                           "all")
                                                           
                except TypeError:
                    continue
                try:
                    dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                except KeyError:
                    continue
                length = utils.get_length(my_file)
                try:
                    dist, branch, delay = utils.fit_distance(conc_dict["Ca"],
                                                             dt,
                                                             method="regular",
                                                             length=length)
                                                                    
                except TypeError:
                    continue
                out_data[k, j, i] = delay.mean()
                out_data_error[k, j, i] = delay.std()/len(delay)**0.5


        my_frame = pd.DataFrame(out_data[:,j,:], index=types,
                                columns=stims)
        my_frame_error = pd.DataFrame(out_data_error[:,j,:],
                                      index=types, columns=stims)
        my_frame.to_csv(save_fname % ("data", diam))
        my_frame_error.to_csv(save_fname % ("data_error", diam))
        print(diam)
        ctrl = out_data[0, j, :]
        for k, (d, fname) in enumerate(directories[1:]):
            print(os.path.join(d, fname))
            paradigm = out_data[k+1, j, :]
            error = out_data_error[k+1, j, :]
            print(ctrl-paradigm)
            print(sum((paradigm - ctrl)/error)/sum(error))
            
        
