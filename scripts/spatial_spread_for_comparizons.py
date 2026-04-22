import sys
import os

import numpy as np

import utility_functions as utils


colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}
directories = [
    [
        "Ca_wave_RyR2CaM_simple_SERCA_SOCE_%s",
        "model_RyR2CaM%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_simple_SERCA_SOCE_%s",
        "model_RyR%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    [
        "Ca_wave_no_RyR_simple_SERCA_SOCE_%s",
        "model_noRyR%s_simple_SERCA_SOCE_tubes_diam_%s_um_50_um_%s_nM.h5"
    ],
    ["Ca_wave_normal_SERCA_aging_%s",
     "model_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
    ["Ca_wave_RyR2CaM_aging_%s",
     "model_RyR2CaM_aging%s_simple_SERCA_tubes_diam_%s_um_50_um_%s_nM.h5"],
]
stims = ["0175", "0350", "0700", "1050", "2000"]
dend_diam = ["1.2", "2.4", "6.0"]
output_name = "all"
if __name__ == "__main__":
    spatial_spread_fname = "spatial_extent_4_cond.csv"
    with open(spatial_spread_fname, "w") as f:
        f.write("condition,model,stim,y\n")
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
                    x_err = []
                    for i, stim in enumerate(stims):
                        #print(fname)
                        new_fname = fname % (stim_type, diam, stim)
                        my_file = os.path.join(my_path % diam, new_fname)
                        try:
                            conc_dict, times_dict = utils.get_conc(my_file,
                                                                   ["Ca"],
                                                                   reg_list,
                                                                   output_name)
                        except TypeError:
                            continue
                        dt = times_dict["trial0"][1]-times_dict["trial0"][0]
                        length = utils.get_length(my_file)
                        dist, branch, delay = utils.fit_distance(conc_dict["Ca"],
                                                                 dt,
                                                                 length=length,
                                                                 find_middle=True)
                        for l, b in enumerate(branch):
                            f.write("%s,%s,%5.3f,%5.3f\n" %(d[:-3], diam,
                                                            b/1000, delay[l]))
                            # print("%s,%s,%s,%4.2f" %(d[:-3], diam,
                            #                              stim, delay[l]))
