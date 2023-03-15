import subprocess
import numpy as np


model_text = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<SDRun xmlns:xi="http://www.w3.org/2001/XInclude" xmlns="http://stochdiff.textensor.org">
    <xi:include href="../Rxn_module_RyR.xml" />
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
    <runtime>1000</runtime>

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



IC_text = """<?xml version="1.0" encoding="utf-8"?>
<InitialConditions>
  <ConcentrationSet>
    <NanoMolarity specieID="Ca" value="%f"/>
    <NanoMolarity specieID="RyRC1"      value="0.1"    />
  </ConcentrationSet>
</InitialConditions>
"""

ca_conc_file = "po_pCa.csv"
exp_res = np.loadtxt(ca_conc_file, skiprows=1, delimiter=',')
ca_conc_list = exp_res[:, 0]

for ca_conc in ca_conc_list:
    ca_conc_nM = int(np.ceil(ca_conc*1e9))
    IC_name = "Ca_%d.xml" % ca_conc_nM
    model_name = "RyR_model_Ca_%d.xml" % ca_conc_nM
    fic = open(IC_name, "w")
    fic.write(IC_text % ca_conc_nM)
    fic.close()
    fm = open(model_name, "w")
    fm.write(model_text % IC_name)
    fm.close()
    
