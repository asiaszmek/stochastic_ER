import os
import subprocess


file_list = os.listdir()
string1 = "Stim_middle_40_ms"
string2 = "Stim_middle_3_ms"

for fname in file_list:
    if fname.startswith(string1) and fname.endswith(".xml"):
        split = fname.split(string1)
        new_fname = string2 + split[-1]
        subprocess.run(["cp", fname, new_fname], capture_output=True)
