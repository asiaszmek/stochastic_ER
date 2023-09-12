import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "700" in fname and fname.endswith(".xml"):
        split = fname.split("700")
        new_fname = split[0]+"0700"+split[-1]
        print(new_fname)
        subprocess.run(["git", "mv", fname, new_fname], capture_output=True)
