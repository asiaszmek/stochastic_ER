import os
import sys
import h5py
from skimage.transform import hough_line, hough_line_peaks, probabilistic_hough_line

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import utility_functions as ut
colors =  {
    "1.2": 'tab:blue',
    "2.4": 'tab:purple',
    "6.0": 'tab:green'
}

names_dict = {
    # "ctrl" : "model_RyRCaM_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    # "normal PMCA + no RyR2": "model_noRyR_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    # "normal PMCA + RyR2": "model_RyR_simple_SERCA_SOCE_tubes_diam_%s_um_2_um_dendrite.h5",
    "low PMCA + RyR2CaM": "model_RyRCaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
   "low PMCA, RyR2 no CaM": "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "low PMCA + 50% RyR2 + 50% RyR2CaM": "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
   "low PMCA + 2x(50% RyR2 + 50% RyR2CaM)": "model_2x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5"
}

dend_diam = ["1.2", "2.4", "6.0"]
output = "__main__"

def adjust_axes(ax):
    mini = min([min(x.get_xlim()) for x in ax])
    maxi = max([max(x.get_xlim()) for x in ax])
    for x in ax:
      
        x.set_xlim([mini, maxi])

        
def find_coords(P0, P1):
    return [[min(P0[0], P1[0]), min(P0[1], P1[1])],
            [max(P0[0], P1[0]), max(P0[1], P1[1])]]



def is_y(miny, maxy, old_miny, old_maxy):
    if miny <= old_maxy and miny >= old_miny :
        return True
    if maxy <= old_maxy and maxy >= old_miny:
        return True
    

if __name__ == "__main__":
    data_dir = os.path.join("..", "stacked_ER")
    x_labels = list(names_dict.keys())
    # fig_m_ca, ax_m_ca = plt.subplots(1, len(dend_diam),
    #                                   figsize=(len(dend_diam)*5, 5))
    # fig_peak_no, ax_peak_no = plt.subplots(1, len(dend_diam),
    #                                        figsize=(len(dend_diam)*5, 5))
    # fig_m_peak_amp, ax_m_peak_amp = plt.subplots(1, len(dend_diam),
    #                                              figsize=(len(dend_diam)*5, 5))
    # fig_m_peak_width, ax_m_peak_width = plt.subplots(1, len(dend_diam),
    #                                                  figsize=(len(dend_diam)*5,
    #                                                           5)) 
    # fig_m_peak_dist, ax_m_peak_dist = plt.subplots(1, len(dend_diam),
    #                                                figsize=(len(dend_diam)*5,
    #                                                         5))
   
    
    for i, d in enumerate(dend_diam):
        no_mean = []
        no_std = []
        amp_mean = []
        amp_std = []
        width_mean = []
        width_std = []
        dist_mean = []
        dist_std = []
        for fname in names_dict.values():
            path = os.path.join(data_dir, fname % d)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            peak_no = []
            peak_amp = []
            peak_width = []
            peak_len = []
            peak_dist = []
            peak_start = []
            for trial in ["trial0", "trial1", "trial2", "trial3", "trial4"]:
                try:
                    data = ut.get_populations(my_file, trial=trial,
                                              output=output)
                except KeyError:
                    continue
                specie_idx = ut.get_all_species(my_file,
                                                output=output).index("Ca")

                volume = np.array([g[12] for g in grid_list])
                conc =  ut.nano_molarity(data[:, :, specie_idx],
                                         volume)
                new_conc = np.reshape(conc, (data.shape[0],
                                             102,
                                             len(volume)//102)).mean(axis=2).T
              
                fig, ax = plt.subplots(1, 3, figsize=(15, 5))
                where = np.where(new_conc > 2.5*76)
                if not len(where[0]):
                    continue
                edges = np.zeros_like(new_conc)
                for x, y in zip(where[0], where[1]):
                   edges[x, y] = 1

                

           
                clusters = [[[where[0][0], where[1][0]],
                             [where[0][0], where[1][0]]]]
                indices = zip(where[0][1:], where[1][1:])
                for i in range(len(where[0][1:])):
                    new_x, new_y = where[0][i+1], where[1][i+1]
                    for  j, c in enumerate(clusters):
                        if c[0][0]-10 <= new_x <= c[1][0]+10:
                            if c[0][1]-4 <= new_y <= c[1][1]+4:
                                new_minx = min([c[0][0], c[1][0], new_x])
                                new_maxx = max([c[0][0], c[1][0], new_x])
                                new_miny = min([c[0][1], c[1][1], new_y])
                                new_maxy = max([c[0][1], c[1][1], new_y])
                                clusters[j] = [[new_minx, new_miny],
                                               [new_maxx, new_maxy]]
                                break
                    else:
                        clusters.append([[new_x, new_y], [new_x, new_y]])

                new_clusters = []
             
                for i, c in enumerate(clusters):
                    new_p1, new_p2 = c
                    for j, nc in enumerate(new_clusters):
                        p1, p2 = nc
                        if (p1[0] <= new_p1[0] <= p2[0]
                            or p1[0] <= new_p2[0] <= p2[0]):
                            if (p1[1] <= new_p1[1] <= p2[1]
                                or p1[1]<= new_p2[1] <=p2[1]):
                                c_minx = min([p1[0], p2[0],
                                              new_p1[0], new_p2[0]])
                                c_maxx = max([p1[0], p2[0],
                                              new_p1[0], new_p2[0]])
                                c_miny = min([p1[1], p2[1],
                                              new_p1[1], new_p2[1]])
                                c_maxy = max([p1[1], p2[1],
                                              new_p1[1], new_p2[1]])
                                new_clusters[j] = [[c_minx, c_miny],
                                                   [c_maxx, c_maxy]]
                                break
                    else:   
                        if  c[1][1] - c[0][1]>=5 and c[1][0] - c[0][0]>=5:
                            new_clusters.append(c)
                if len(new_clusters):

                    new_clusters = sorted(new_clusters, key=lambda x:x[0][0])
                    peak_no.append(len(new_clusters))
                    
                    for j, nc in enumerate(new_clusters):
                        p1, p2 = nc
                        peak_width.append((p2[1] - p1[1])*0.2)
                        peak_len.append((p2[0]-p1[0])/2+0.5)
                        peak_amp.append(new_conc[p1[0]:p2[0]+1,
                                             p1[1]:p2[1]+1].mean())
                        
                        peak_start.append(new_conc[:, p1[1]].argmax()/2)
                        if j:
                            o_p1, o_p2 = new_clusters[j-1]
                            new_t = (p2[1]+p1[1])/2*0.2
                            old_t = (o_p2[1]+o_p1[1])/2*0.2
                            peak_dist.append(abs(new_t - old_t))
       
        print(fname)
        peak_no = np.array(peak_no)
        peak_amp = np.array(peak_amp)
        peak_width = np.array(peak_width)
        peak_len = np.array(peak_len)
        peak_dist = np.array(peak_dist)
        peak_start = np.array(peak_start)
        print(peak_no, peak_no.mean(), peak_no.std()/len(peak_no)**.5)
        print(peak_len, peak_len.mean(), peak_len.std()/len(peak_len)**.5)
        print(peak_width, peak_width.mean(), peak_width.std()/len(peak_width)**.5)
        print(peak_dist, peak_dist.mean(), peak_dist.std()/len(peak_dist)**.5)
        print(peak_start, peak_start.mean(), peak_start.std()/len(peak_start)**.5)

        print(peak_amp, peak_amp.mean(), peak_amp.std()/len(peak_amp)**.5)
        
        
