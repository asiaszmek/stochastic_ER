import os
import subprocess


file_list = os.listdir()


for fname in file_list:
    if fname.startswith("model_RyR") and fname.endswith(".xml"):
        split = fname[:-4]
        new_fname = split+"Fura2.xml"
        print(new_fname)
        subprocess.run(["cp", fname, new_fname], capture_output=True)
