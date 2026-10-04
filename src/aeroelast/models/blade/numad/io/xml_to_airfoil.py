import re

import numpy as np


def xml_to_airfoil(airfoil, filecontents: list[str] | None = None):
    """TODO docstring

    Parameters
    ----------

    Returns
    -------
    """
    if filecontents is None:
        raise ValueError("xml_to_airfoil requires the file contents")
    # DEVNOTE: this works for the kinds of af files
    #          found in the BAR reference blades
    #          but not for general xmls
    #          this should be extended to be more general at some point
    # [coords, reference] = readAirfoilXML(filecontents)

    # The following regular expression pattern matches any number of characters
    # found between the opening and closing "reference" tags
    coords = np.empty((0, 2))

    fulltext = "".join(filecontents)

    pattern = "<reference>(.*)</reference>"
    t = re.search(pattern, fulltext)
    if t is None:
        raise ValueError("airfoil XML has no <reference> tag")
    reference = t.group(1)

    for line in filecontents:
        # check if there is a tag
        if re.search("<", line):
            continue
        # otherwise, assume coordinate data
        else:
            if re.search("\t", line):
                line = line.replace("\t", " ")
            x, y = line.split(" ")
            x = float(x)
            y = float(y)
            coords = np.append(coords, [[x, y]], axis=0)
    airfoil.reference = reference
    airfoil.coordinates = coords
