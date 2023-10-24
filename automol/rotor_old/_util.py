"""
 Various utility functions
"""


def sort_tors_names(tors_names):
    """ sort torsional names so that Dn where n is ascending order
    """
    tors_names = list(tors_names)
    tors_names.sort(key=lambda x: int(x.split('D')[1]))
    return tors_names
