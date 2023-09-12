import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "350" in fname and fname.endswith(".xml"):
        split = fname.split("350")
        new_fname = split[0]+"0350"+split[-1]
        print(new_fname)
        subprocess.run(["git", "mv", fname, new_fname], capture_output=True)
