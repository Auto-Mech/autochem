"""
Functions to read data files
"""

import os
import numpy


PATH = os.path.dirname(os.path.realpath(__file__))


def read_file(path_lst, file_name):
    """ Read a file
    """
    file_path = os.path.join(PATH, *path_lst, file_name)
    with open(file_path, 'r') as fobj:
        file_str = fobj.read()
    return file_str


def load_numpy_string_file(path_lst, file_name):
    """ Read a file with numpy
    """
    file_path = os.path.join(PATH, *path_lst, file_name)
    file_lst = numpy.loadtxt(file_path, dtype=str)

    return file_lst
