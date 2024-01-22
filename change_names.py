import os
import glob
import subprocess


directories = glob.glob("Ca_wave*")
print(directories)
for directory in directories:
    file_list = os.listdir(directory)
    

    for fname in file_list:
        if "1050" in fname and fname.endswith(".xml"):
            split = fname.split("1050")
            new_fname = os.path.join(directory, split[0]+"2000"+split[-1])
            print(os.path.join(directory, fname), new_fname)
            subprocess.run(["cp","-v", os.path.join(directory, fname), new_fname],
                           capture_output=True)
