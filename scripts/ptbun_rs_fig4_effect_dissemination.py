import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import utility_functions as utils


colors = {"no_SOCE_no_CaM": "tab:blue",
          "_no_CaM": "tab:green",
          "no_SOCE_CaM": "tab:cyan",
          "_CaM": "tab:olive",
}
labels = {"no_SOCE_no_CaM": "no SOCE RyR dis-inh.",
          "_no_CaM": " RyR dis-inh.",
          "no_SOCE_CaM": "no SOCE",
          "_CaM": "ctrl",
}
base = "dend"
reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
            "dend06", "dend07", "dend08", "dend09",]
for i in range(10, 102, 1):
    reg_list.append("%s%d" %(base, i))

file_path = os.path.abspath(__file__)
list_fp = os.path.split(file_path)
cur_dir = os.path.join(os.path.join(list_fp[0], ".."))
basic_RyR_no_SOCE_dir = "Ca_wave_simple_SERCA_no_SOCE_Breit_2018"
basic_RyR_SOCE_dir = "Ca_wave_simple_SERCA_SOCE"
RyRCaM_no_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE"
RyRCaM_SOCE_dir = "Ca_wave_RyR2CaM_simple_SERCA_SOCE"
output_name = "all"
stim = "0700"
stim_label = "4 uM Ca injection"
branch_diam = 1.2
t_start = 3000
idx_start = t_start
base_noSOCE = "model_RyR_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
base_SOCE = "model_RyR_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

CaM_noSOCE = "model_RyR2CaM_simple_SERCA_tubes_diam_%2.1f_um_50_um_%s_nM.h5"
CaM_SOCE = "model_RyR2CaM_simple_SERCA_SOCE_tubes_diam_%2.1f_um_50_um_%s_nM.h5"

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

mini = []
maxi = []


fname = base_noSOCE % (branch_diam, stim)
full_name_noSOCE_noCaM = os.path.join(cur_dir, basic_RyR_no_SOCE_dir, fname)
voxels, time, ca_noSOCE_noCaM = utils.get_conc(full_name_noSOCE_noCaM, ["Ca"], reg_list, output_name)

fname_SOCE = base_SOCE % (branch_diam, stim)
full_name = os.path.join(cur_dir, basic_RyR_SOCE_dir, fname_SOCE)
voxels, time, ca_SOCE_noCaM = utils.get_conc(full_name, ["Ca"], reg_list, output_name)

fname = CaM_noSOCE % (branch_diam, stim)
full_name = os.path.join(cur_dir, RyRCaM_no_SOCE_dir, fname)
voxels, time, ca_noSOCE_CaM = utils.get_conc(full_name, ["Ca"], reg_list, output_name)

        
fname_SOCE = CaM_SOCE % (branch_diam, stim)
full_name = os.path.join(cur_dir, RyRCaM_SOCE_dir, fname_SOCE)
voxels, time, ca_SOCE_CaM = utils.get_conc(full_name, ["Ca"], reg_list, output_name)
v_raw_min = min(ca_SOCE_CaM.min(),  ca_SOCE_noCaM.min())
v_raw_max = max(ca_SOCE_CaM.max(),  ca_SOCE_noCaM.max())
# SOCE no CaM -- basal
image_full = ax[0][0].imshow(ca_SOCE_CaM.T, aspect="auto",
                             interpolation="none",
                             origin="lower", extent = [time[0]*1e-3,
                                                       time[-1]*1e-3,
                                                       voxels[0],
                                                       voxels[-1]],
                             vmin=v_raw_min,
                             vmax=v_raw_max,
                             cmap=plt.get_cmap("Reds"))
fig.colorbar(image_full, ax=ax[0][0])
ax[0][0].set_title("control")
ax[0][0].set_ylabel("x [um]")

soce_alone =  ca_SOCE_noCaM - ca_noSOCE_noCaM
mini.append(soce_alone.min())
maxi.append(soce_alone.max())
cam_alone = ca_SOCE_noCaM - ca_SOCE_CaM
mini.append(cam_alone.min())
maxi.append(cam_alone.max())

vmin = min(mini)
vmax = max(maxi)
if abs(vmin) > abs(vmax):
    vmax = abs(vmin)
else:
    vmin = -vmax


image_soce = ax[0][1].imshow(soce_alone.T, aspect="auto",
                             interpolation="none",
                             origin="lower", extent = [time[0]*1e-3,
                                                       time[-1]*1e-3,
                                                       voxels[0],
                                                       voxels[-1]],
                             vmin=vmin, vmax=vmax,
                             cmap=plt.get_cmap("bwr"))
ax[0][1].set_title("SOCE contribution")
fig.colorbar(image_soce, ax=ax[0][1])
image_soce = ax[1][0].imshow(cam_alone.T, aspect="auto",
                             interpolation="none",
                             origin="lower", extent = [time[0]*1e-3,
                                                       time[-1]*1e-3,
                                                       voxels[0],
                                                       voxels[-1]],
                             vmin=vmin, vmax=vmax,
                             cmap=plt.get_cmap("bwr"))
ax[1][0].set_title("RyR dis-inhibition contribution")
fig.colorbar(image_soce, ax=ax[1][0])
image_disease = ax[1][1].imshow(ca_noSOCE_noCaM.T, aspect="auto",
                                interpolation="none",
                                origin="lower", extent = [time[0]*1e-3,
                                                          time[-1]*1e-3,
                                                          voxels[0],
                                                          voxels[-1]],
                                vmin=v_raw_min,
                                vmax=v_raw_max,
                                
                                cmap=plt.get_cmap("Reds"))
ax[1][1].set_title("RyR dis-inhibition")
print(ca_SOCE_noCaM.T.max(), ca_SOCE_CaM.T.max())
fig.colorbar(image_disease)
ax[1][0].set_ylabel("x [um]")
ax[1][0].set_xlabel("time [sec]")
ax[1][1].set_xlabel("time [sec]")


fig.savefig("Ca_mechanisms_contribution.svg", dpi=100, bbox_inches="tight")
plt.show()
