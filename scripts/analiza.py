#!/usr/bin/env python
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import utility_functions as utils


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    for fname in sys.argv[1:]:
        my_file = h5py.File(fname, 'r')
        output_dict = utils.get_output_regions(my_file)
        for trial in my_file.keys():
            if trial != "model":
                utils.save_concentrations(my_file, fname[:-3], '__main__',
                                    trial=trial)

    print('Done')

