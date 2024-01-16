import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "largerER_largerER" in fname and fname.endswith(".xml"):
        split = fname.split("largerER_largerER")
        new_fname = split[0]+"largerER"+split[-1]
        print(new_fname)
        subprocess.run(["mv", "-v", fname, new_fname], capture_output=True)
