import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if "Fura2" in fname and fname.endswith(".xml"):
        split = fname.split(".xmlFura2")
        new_fname = split[0]+"Fura2"+split[-1]
        print(new_fname)
        subprocess.run(["git", "mv",  fname, new_fname], capture_output=True)
