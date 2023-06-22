import sys
import subprocess
from datetime import date
from lxml import etree
import h5py
import numpy as np
import matplotlib.pyplot as plt

ca_conc_file = "xu_meissner_RyR2CaM_po.csv"
ryr_op_fname = "xu_meissner_RyR2CaM_minopentime.csv"
ryr_cl_fname = "xu_meissner_RyR2CaM_minclosedtime.csv"

model_text = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<SDRun xmlns:xi="http://www.w3.org/2001/XInclude" xmlns="http://stochdiff.textensor.org">
    <xi:include href="%s" />
    <xi:include href="Morph.xml" />
    <xi:include href="%s" />
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
    <runtime>1000000</runtime>

    <!-- set the seed to get the same spines each time testing -->
    <spineSeed>123</spineSeed>

    <discretization>
      <defaultMaxElementSide>10</defaultMaxElementSide>
      <surfaceLayers>10</surfaceLayers> 
    </discretization>
    <tolerance>0.01</tolerance>

    <outputInterval>100</outputInterval>

    <calculation>GRID_ADAPTIVE</calculation>

</SDRun>"""



IC_text = """<?xml version="1.0" encoding="utf-8"?>
<InitialConditions>
  <ConcentrationSet>
    <NanoMolarity specieID="Ca" value="%f"/>
    <NanoMolarity specieID="CaM" value="1000"/>
    <NanoMolarity specieID="RyR"      value="10"    />
  </ConcentrationSet>
</InitialConditions>
"""


Rxn_file = "Rxn.xml"


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

     

if __name__ == "__main__":

        
    for i, ca_conc in enumerate([0, 10, 50, 100, 500, 1000, 5000, 10000, 20000, 100000]):
        ca_conc_nM = 2*int(np.ceil(ca_conc))
        IC_name = "Ca_%d.xml" % ca_conc_nM
        model_name = "RyR2_model_Ca_%d.xml" % ca_conc_nM
        output_name = "RyR2_model_Ca_%d.h5" % ca_conc_nM
        fic = open(IC_name, "w")
        fic.write(IC_text % ca_conc_nM)
        fic.close()
        fm = open(model_name, "w")
        fm.write(model_text % (Rxn_file, IC_name))
        fm.close()
        process = subprocess.run(["/usr/lib/jvm/java-8-openjdk-amd64/bin/java",
                                  "-jar",
                                  "/home/jszmek/new_neurord/neurord-3.2.3-all-deps.jar",
                                  "-Dneurord.trials=4", model_name],
                                 capture_output=True)
        print(process.returncode)
