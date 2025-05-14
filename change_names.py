import os
import glob
import subprocess


directory = "Ca_wave_RyR2CaM_simple_SERCA_no_SOCE_largerER"
file_list = glob.glob(os.path.join(directory, "*xml"))
print(file_list)    

for fname in file_list:
  
    if "model" in fname:
        
        split = fname.split("_no_SOCE_largerER")
        new_name = os.path.join("Ca_wave_RyR2CaM_simple_SERCA_SOCE_largerER",
                                split[-1][1:])
        
        print(os.path.join(directory, fname), new_name)
        subprocess.run(["git","mv", fname,
                        new_name],
                       capture_output=True)
