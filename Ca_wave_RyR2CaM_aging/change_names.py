import os
import glob
import subprocess



file_list = os.listdir(".")
    

for fname in file_list:
    if "RyR2CaM" not in fname and fname.endswith(".xml"):
        split = fname.split("aging")
        new_fname = split[0]+"RyR2CaM_aging"+split[-1]
        print(fname, new_fname)
        subprocess.run(["cp","-v", fname, new_fname],
                       capture_output=True)
