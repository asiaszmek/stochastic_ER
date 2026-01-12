import os
import glob
import subprocess

if __name__ == "__main__":

    lista = glob.glob("Ca_wave*")

    for item in lista:
        zip_dest = item + ".zip"
        subprocess.run(["zip", "-r", zip_dest, item],
                       capture_output=True)
