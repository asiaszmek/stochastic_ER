import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "Fura2" in fname and fname.endswith(".h5"):
        split = fname.split("_SERCA_SOCE")
        
        new_fname = split[0]+"_SERCA_SOCE_Fura2"+split[-1][:-8]+".h5"
        print(new_fname)
        subprocess.run(["mv", fname, new_fname], capture_output=True)
