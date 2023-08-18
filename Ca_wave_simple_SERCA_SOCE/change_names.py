import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "1.2" in fname:
        continue
    if fname.startswith("model_RyR") and fname.endswith(".xml") and "tubes" in fname:
        split = fname.split("tubes")
        new_fname = split[0]+"nc_tubes"+split[1]
        print(new_fname)
        subprocess.run(["cp", fname, new_fname], capture_output=True)
