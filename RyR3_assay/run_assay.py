import sys
import subprocess
from lxml import etree
import h5py
import numpy as np
import matplotlib.pyplot as plt

KL_list = ["RyR", "RyRC1",
           "RyRC2", "RyRC3",
           "RyRC4", "RyRC5",
           "RyRI"]

def get_key(cell):
    if cell[18]:
        return cell[15].decode('utf-8') + '_' + cell[18].decode('utf-8')
    return cell[15].decode('utf-8')

def get_regions(my_file):
    grid_list = get_grid_list(my_file)
    return sorted(list(set([get_key(grid) for grid in grid_list])))

def region_volumes(my_file):
    grid_list = get_grid_list(my_file)
    regions = get_regions(my_file)
    volumes = {}
    for region in regions:
        volumes[region] = 0
    for cell in grid_list:
        key = get_key(cell)
        volumes[key] += float(cell[12])
    return volumes

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

def get_output_regions(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    outputs = {}
    for son in root:
        if son.tag.endswith('OutputScheme'):
            for grandson in son:
                outputs[grandson.get("filename")] = grandson.get("region")
    return outputs

def sum_volume(my_file, region_list):
    grid_list = get_grid_list(my_file)
    vol_sum = 0
    volumes = region_volumes(my_file)
    for region in region_list:
        if region in volumes:
            vol_sum += volumes[region]
    return vol_sum

def get_all_closed(data, species, specie_list):
    count = 0
    sum_times = 0
    out = []
    for specie in specie_list:
        specie_state = data[:, 0, species.index(specie)]
        if "RyR" not in specie_list:
            count += len(np.where((specie_state[1:] - specie_state[0:-1])==1)[0])
    
        sum_times += specie_state.sum()
        out.append(specie_state[-1])
    return sum_times, count, out

def get_numbers(my_file, output="all"):
    output_dict = get_output_regions(my_file)
    Ca_conc = []
    open_ryr3 = []
    mean_c_t = []
    mean_o_t = []
    for trial in my_file.keys():
        if trial == "model":
            continue
        times = get_times(my_file, trial=trial, output=output)
        vol = sum_volume(my_file, ["dend"])
        species = get_all_species(my_file, output=output)
        data = get_populations(my_file, trial=trial, output=output)
        dt = times[1]-times[0]
        exp_len = int((times[-1])/dt)
        mean_ca = data[:, 0, species.index("Ca")].mean()*10/6.023/vol
        ryr_basal = data[0, 0, species.index("RyRC1")]
        open_sum, tot_no, ends = get_all_closed(data, species, ["RyRO1", "RyRO2"])
        if 1 in ends:
            end_closed = False
        else:
            end_closed = True
        p_open_ryr = open_sum/ryr_basal/exp_len
        Ca_conc.append(mean_ca)
        open_ryr3.append(p_open_ryr)

        if ryr_basal != 1:
            continue
        if tot_no > 0:
            mean_o_t.append(dt*open_sum/tot_no)     
        
        sum_closed, tot_nc, ends  = get_all_closed(data, species, ["RyRC1", "RyRC2"] )

        if end_closed:
            tot_nc += 1
        tot_nc += tot_no
        if tot_nc > 0:
            mean_c_t.append(dt*sum_closed/tot_nc)
    return Ca_conc, open_ryr3, mean_o_t, mean_c_t



        
model_text = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<SDRun xmlns:xi="http://www.w3.org/2001/XInclude" xmlns="http://stochdiff.textensor.org">
    <xi:include href="../Rxn_module_RyR_KeizerLevine.xml" />
    <xi:include href="Morph.xml" />
    <xi:include href="%s" />
    <xi:include href="IO_Ca.xml"/>
    <!--2D means the morphology is interpreted like a flatworm, 3D for
roundworms. The 2D case is good for testing as it is easy to visualize the
results (also, 3D may not work yet...)  -->
   
    <geometry>          2D           </geometry>
    <depth2D>           10          </depth2D>
    <distribution>      BINOMIAL     </distribution>
    <algorithm>         INDEPENDENT  </algorithm>
    <simulationSeed>    245         </simulationSeed>
    <outputQuantity>NUMBER</outputQuantity>

    <!-- run time for the calculation, milliseconds -->
    <runtime>143000</runtime>

    <!-- set the seed to get the same spines each time testing -->
    <spineSeed>123</spineSeed>

    <discretization>
      <defaultMaxElementSide>10</defaultMaxElementSide>
      <surfaceLayers>10</surfaceLayers> 
 
      

    </discretization>


    <!-- the tolerace is not used yet -->
    <tolerance>0.001</tolerance>

    <outputInterval>500</outputInterval>



    <calculation>GRID_ADAPTIVE</calculation>

</SDRun>"""



IC_text_KL = """<?xml version="1.0" encoding="utf-8"?>
<InitialConditions>
  <ConcentrationSet>
    <NanoMolarity specieID="Ca" value="%f"/>
    <NanoMolarity specieID="RyRC1"      value="0.1"    />
  </ConcentrationSet>
</InitialConditions>
"""

IO_KL = """
<OutputScheme>
  <OutputSet filename = "all"  dt=".1">
    <OutputSpecie name="Ca"/>
    <OutputSpecie name="RyRO1"/>
    <OutputSpecie name="RyRO2"/>
    <OutputSpecie name="RyRC1"/>
    <OutputSpecie name="RyRC2"/>
  </OutputSet>

</OutputScheme>"""


IO_Dura = """
<OutputScheme>
  <OutputSet filename = "all"  dt=".1">
    <OutputSpecie name="Ca"/>
    <OutputSpecie name="RyRO1"/>
    <OutputSpecie name="RyRO2"/>
    <OutputSpecie name="RyRC1"/>
    <OutputSpecie name="RyRC2"/>
    <OutputSpecie name="RyRC4"/>
    <OutputSpecie name="RyRC5"/>
    <OutputSpecie name="RyR"/>
    <OutputSpecie name="RyRI"/>

  </OutputSet>

</OutputScheme>"""


f = open("IO_Ca.xml", "w")
f.write(IO_KL)
f.close()
ca_conc_file = "po_pCa.csv"
ryr_op_fname = "open_duration_pCa.csv"
ryr_cl_fname = "closed_duration_pCa.csv"
exp_res = np.loadtxt(ca_conc_file, skiprows=1, delimiter=',')
ca_conc_list = exp_res[:, 0]
output = np.zeros(exp_res.shape)
mean_times = []
for i, ca_conc in enumerate(ca_conc_list):

    ca_conc_nM = int(np.ceil(ca_conc*1e9))
    IC_name = "Ca_%d.xml" % ca_conc_nM
    model_name = "RyR_model_Ca_%d.xml" % ca_conc_nM
    output_name = "RyR_model_Ca_%d.h5" % ca_conc_nM
    fic = open(IC_name, "w")
    fic.write(IC_text_KL % ca_conc_nM)
    fic.close()
    fm = open(model_name, "w")
    fm.write(model_text % IC_name)
    fm.close()
    process = subprocess.run(["/usr/lib/jvm/java-8-openjdk-amd64/bin/java",
                              "-jar",
                              "/home/jszmek/new_neurord/neurord-3.2.3-all-deps.jar",
                              "-Dneurord.trials=10", model_name],
                             capture_output=True)
    if not process.returncode:
        my_file = h5py.File(output_name, 'r')
        conc, po, t_o, t_c = get_numbers(my_file, output="all")
        output[i, 0] = np.mean(conc)
        output[i, 1] = np.mean(po)
        if len(t_o) == 0 or len(t_c) == 0:
            continue
        mean_times.append([np.mean(conc), np.mean(t_o), np.mean(t_c)])

            
exp_open = np.loadtxt(ryr_op_fname, skiprows=1, delimiter=",")
exp_closed = np.loadtxt(ryr_cl_fname, skiprows=1, delimiter=",")
mean_times_a = np.array(mean_times)

fig, ax = plt.subplots(1)
ax.set_xscale('log')
ax.plot(exp_res[:, 0], exp_res[:, 1], "d", color="tab:blue", label="experimental data")
ax.plot(output[:, 0]*1e-9, output[:, 1], "d", color="tab:red", label="model data")
ax.legend()
ax.set_xlabel("Concentration [M]")
ax.set_ylabel("RyR3 open probability")

fig, ax = plt.subplots(1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(exp_open[:, 0], exp_open[:, 1], "d", color="tab:blue",
        label="exp open")
ax.plot(exp_closed[:, 0], exp_closed[:, 1], "d", color="tab:green",
        label="exp closed")
ax.plot(mean_times_a[:, 0]*1e-9, mean_times_a[:, 1], "d",
        label="model open", color="tab:cyan")
ax.plot(mean_times_a[:, 0]*1e-9, mean_times_a[:, 2], "d",
        label="exp closed", color="tab:olive")

ax.legend()
ax.set_xlabel("Concentration [M]")
ax.set_ylabel("Time [ms]")


plt.show()
