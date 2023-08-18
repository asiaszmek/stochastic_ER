import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "SOCE_SOCE" in fname and fname.endswith(".xml"):
        split = fname.split("SOCE_SOCE")
        new_fname = split[0]+"SOCE"+split[-1]
        print(new_fname)
        subprocess.run(["git", "mv", fname, new_fname], capture_output=True)
