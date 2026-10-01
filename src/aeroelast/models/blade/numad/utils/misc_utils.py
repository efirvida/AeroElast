########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import numpy as np


def setup_logging(file_name):
    import logging

    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    fh = logging.FileHandler(file_name + ".log", mode="w")
    log.addHandler(fh)
    return log


def full_keys_from_substrings(component_dict, key_list):
    return [
        key
        for key in component_dict
        if all(substring.lower() in key.lower() for substring in key_list)
    ]


# SED-like substitution


def _parse_data(data):
    """Helper function for parsing data from blade yaml files.

    Parameters
    ----------
    data
        a number or list of numbers where numbers can be floats or strings
        e.g. 3.0 or '1.2e2'

    Returns
    -------
    parsed_data
        single float or array of floats
    """
    try:
        # detect whether data is list
        _ = data + []
    except TypeError:  # case for single data point
        parsed_data = float(data)
    else:
        parsed_data = np.array([float(val) for val in data])  # case for list of data points
    return parsed_data
