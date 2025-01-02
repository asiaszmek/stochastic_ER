import os
import sys
import h5py
import numpy as np
from scipy.fft import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
import utility_functions as ut
colors =  {
    "1.2": 'tab:blue',
    "2.4": 'tab:purple',
    "6.0": 'tab:green'
}

names_dict = {
    "ctrl" :
    "model_RyRCaM_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    "normal PMCA + no RyR2":
    "model_noRyR_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    "normal PMCA + RyR2":
    "model_RyR_simple_SERCA_SOCE_tubes_diam_%s_um_2_um_dendrite.h5",
    "low PMCA + RyR2CaM":
    "model_RyRCaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "low PMCA, RyR2 no CaM":
    "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "low PMCA + 50% RyR2 + 50% RyR2CaM":
    "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
    "low PMCA + 2x(50% RyR2 + 50% RyR2CaM)":
    "model_2x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5"
}

dend_diam = ["1.2", "2.4", "6.0"]
#trial = "trial0"
output = "__main__"

def adjust_axes(ax):
    mini = min([min(x.get_xlim()) for x in ax])
    maxi = max([max(x.get_xlim()) for x in ax])
    for x in ax:
      
        x.set_xlim([mini, maxi])


def find_peaks(array):
    peak_time = []
    peak_width = []
    peak_amp = []
    if len(np.where(array > 150)[0]):
        peaks = np.where(array > 110)[0]
        if len(peaks):
            ranges = np.where(peaks[1:]-peaks[:-1] > 50)[0]
            idx_pre = peaks[0]
            for r in ranges:
                idx_post = peaks[r]
                idx = np.argmax(array[idx_pre:idx_post])+idx_pre
                peak_amp.append(array[idx])
                peak_time.append(idx*200)
                peak_width.append((idx_post-idx_pre)*200)
                idx_pre = idx_post
            if array[idx_pre:peaks[-1]].max() > 110:
                idx = np.argmax(array[idx_pre:peaks[-1]])+idx_pre
                peak_time.append(idx*200)
                peak_width.append((peaks[-1]-idx_pre)*200)
                peak_amp.append(array[idx])
    return peak_time, peak_width, peak_amp

if __name__ == "__main__":
    data_dir = os.path.join("..", "stacked_ER")
    x_labels = []
    x_labels_pmca = []
    fig_m_ca, ax_m_ca = plt.subplots(1, len(dend_diam),
                                           figsize=(len(dend_diam)*5, 5))
    fig_peak_no, ax_peak_no = plt.subplots(1, len(dend_diam),
                                           figsize=(len(dend_diam)*5, 5))
    fig_m_peak_amp, ax_m_peak_amp = plt.subplots(1, len(dend_diam),
                                                 figsize=(len(dend_diam)*5, 5))
    fig_m_peak_width, ax_m_peak_width = plt.subplots(1, len(dend_diam),
                                                     figsize=(len(dend_diam)*5,
                                                              5)) 
    fig_m_peak_dist, ax_m_peak_dist = plt.subplots(1, len(dend_diam),
                                                   figsize=(len(dend_diam)*5,
                                                            5))
   
    
    for i, d in enumerate(dend_diam):
        means = []
        stds = []
        no_mean = []
        no_std = []
        amp_mean = []
        amp_std = []
        width_mean = []
        width_std = []
        dist_mean = []
        dist_std = []
        x_labels = []
        x_labels_pmca = []
        for key, fname in names_dict.items():
            path = os.path.join(data_dir, fname % d)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            conc_list = []
            conc_std = []
            peak_no = []
            peak_amp = []
            peak_width = []
            peak_dist = []
            for trial in ["trial0", "trial1", "trial2", "trial3"]:
                data = ut.get_populations(my_file, trial=trial, output=output)
                specie_idx = ut.get_all_species(my_file,
                                                output=output).index("Ca")
                volume = sum(list(ut.region_volumes(my_file).values()))
                
                tot_conc =  ut.nano_molarity(data[:, :, specie_idx].sum(axis=1),
                                         volume)
                
                peaks, peak_w, peak_amp_i = find_peaks(tot_conc)

                if len(peaks):
                    peak_no.append(len(peaks))
                    peak_width += peak_w
                    peak_amp += peak_amp_i
                    peak_dist += [d-peaks[i-1] for i, d in enumerate(peaks[1:])]
                else:
                    peak_no.append(0)
                    peak_amp.append(0)
                    peak_width.append(0)
                    peak_dist.append(0)
                conc_list.append(tot_conc.mean())
                conc_std.append(tot_conc.std())

            means.append(np.mean(conc_list))
            stds.append(sum([s**2 for s in conc_std])**0.5/2)
            x_labels.append(key)
            if "low PMCA" in key:
                print(key)
                x_labels_pmca.append(key)
                no_mean.append(np.mean(peak_no))
                no_std.append(np.std(peak_no)/len(peak_no)**0.5)
                amp_mean.append(np.mean(peak_amp))
                amp_std.append(np.std(peak_amp)/len(peak_amp)**.5)
                width_mean.append(np.mean(peak_width)/1000)
                width_std.append(np.std(peak_width)/len(peak_width)**.5/1000)
                dist_mean.append(np.mean(peak_dist)/1000)
                dist_std.append(np.std(peak_dist)/len(peak_dist)**.5/1000)
        print(d, no_mean,  amp_mean, width_mean, dist_mean)
        ax_m_ca[i].errorbar(y=x_labels, x=means, xerr=stds,
                        color=colors[d],
                        marker="o", linestyle="")
        ax_m_ca[i].set_title("diam %s um" % d)
        ax_m_ca[i].set_yticklabels([])
        ax_m_ca[i].set_xlabel("mean Ca (nM)")
        
        ax_peak_no[i].errorbar(y=x_labels_pmca, x=no_mean,
                               xerr=no_std,
                               color=colors[d],
                               marker="o", linestyle="")
        ax_peak_no[i].set_title("diam %s um" % d)
        ax_peak_no[i].set_yticklabels([])
        ax_peak_no[i].set_xlabel("# calcium peaks")
        ax_m_peak_amp[i].errorbar(y=x_labels_pmca, x=amp_mean,
                                  xerr=amp_std,
                                  color=colors[d],
                                  marker="o", linestyle="")
        
        ax_m_peak_amp[i].set_yticklabels([])
        ax_m_peak_amp[i].set_title("diam %s um" % d)
        ax_m_peak_amp[i].set_xlabel("mean Ca peak amplitude (nM)")
        ax_m_peak_width[i].errorbar(y=x_labels_pmca, x=width_mean,
                                    xerr=width_std,
                                    color=colors[d],
                                    marker="o", linestyle="")
        
        ax_m_peak_width[i].set_yticklabels([])
        ax_m_peak_width[i].set_title("diam %s um" % d)
        ax_m_peak_width[i].set_xlabel("mean Ca peak width (s)")

        ax_m_peak_dist[i].errorbar(y=x_labels_pmca, x=dist_mean,
                                   xerr=dist_std,
                                   color=colors[d],
                                   marker="o", linestyle="")
        
        ax_m_peak_dist[i].set_yticklabels([])
        ax_m_peak_dist[i].set_title("diam %s um" % d)
        ax_m_peak_dist[i].set_xlabel("mean interval between peaks (s)")

        
    ax_m_ca[0].set_yticklabels(x_labels)
    ax_peak_no[0].set_yticklabels(x_labels_pmca)
    ax_m_peak_amp[0].set_yticklabels(x_labels_pmca)
    ax_m_peak_width[0].set_yticklabels(x_labels_pmca)
    ax_m_peak_dist[0].set_yticklabels(x_labels_pmca)
    ax_peak_no[0].set_yticklabels(x_labels_pmca)
    ax_m_peak_amp[0].set_yticklabels(x_labels_pmca)
    ax_m_peak_width[0].set_yticklabels(x_labels_pmca)
    ax_m_peak_dist[0].set_yticklabels(x_labels_pmca)
    adjust_axes(ax_m_ca)
    adjust_axes(ax_peak_no)
    adjust_axes(ax_m_peak_amp)
    adjust_axes(ax_m_peak_width)
    adjust_axes(ax_m_peak_dist)
    fig_m_ca.savefig("mean_basal_ca.png", dpi=100,
                 bbox_inches="tight")
    fig_peak_no.savefig("mean_peak_no.png", dpi=100,
                 bbox_inches="tight")
    fig_m_peak_amp.savefig("mean_peak_amp.png", dpi=100,
                 bbox_inches="tight")
    fig_m_peak_dist.savefig("mean_peak_dist.png", dpi=100,
                            bbox_inches="tight")
    fig_m_peak_width.savefig("mean_peak_width.png", dpi=100,
                            bbox_inches="tight")

    plt.show()
