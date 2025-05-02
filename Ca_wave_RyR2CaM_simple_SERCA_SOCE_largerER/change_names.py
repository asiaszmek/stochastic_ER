import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "largerER_sim" in fname and fname.endswith(".xml"):
        split = fname.split("SERCA")
        new_fname = split[0]+"SERCA_SOCE"+split[-1]
        print(new_fname)
        subprocess.run(["git", "mv", "-v", fname, new_fname], capture_output=True)
