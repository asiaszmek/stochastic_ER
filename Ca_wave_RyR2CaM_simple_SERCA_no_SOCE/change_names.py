import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "Fura2" in fname and fname.endswith(".xml"):
        split = fname.split("_SERCA")
        
        new_fname = split[0]+"_SERCA_Fura2"+split[-1][:-9]+".xml"
        print(new_fname)
        subprocess.run(["git", "mv", fname, new_fname], capture_output=True)
