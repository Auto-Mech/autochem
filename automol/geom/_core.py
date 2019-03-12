""" core library defining the geom data structure
"""
from ..constructors.geom import from_data as _from_data


def symbols(geo):
    """ atomic symbols
    """
    if geo:
        syms, _ = zip(*geo)
    else:
        syms = ()
    return syms


def coordinates(geo):
    """ atomic coordinates
    """
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    return xyzs


def set_coordinates(geo, xyz_dct):
    """ set coordinate values for the geometry, using a dictionary by index
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)

    natms = len(syms)
    assert all(idx in range(natms) for idx in xyz_dct)

    xyzs = [xyz_dct[idx] if idx in xyz_dct else xyz
            for idx, xyz in enumerate(xyzs)]
    return _from_data(syms, xyzs)


def is_valid(geo):
    """ is this a valid geometry?
    """
    ret = hasattr(geo, '__iter__')
    if ret:
        ret = all(hasattr(obj, '__len__') and len(obj) == 2 for obj in geo)
        if ret:
            syms, xyzs = zip(*geo)
            try:
                _from_data(syms, xyzs)
            except AssertionError:
                ret = False
    return ret
