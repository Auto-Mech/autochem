""" cartesian geometry writers
"""


def write(syms, xyzs):
    """ write a geometry to a string
    """
    natms = len(syms)
    assert len(xyzs) == natms

    geo_str = '\n'.join('{:2s} {:10.6f} {:10.6f} {:10.6f}'.format(sym, *xyz)
                        for sym, xyz in zip(syms, xyzs))
    return geo_str


def write_xyz(syms, xyzs, comment=None):
    """ write a geometry to a .xyz string
    """
    comment = '' if comment is None else comment

    natms = len(syms)
    assert len(xyzs) == natms

    geo_str = write(syms=syms, xyzs=xyzs)
    xyz_str = ' {:d}\n{:s}\n{:s}'.format(natms, comment, geo_str)
    return xyz_str


def write_xyz_trajectory(syms, xyzs_lst, comments=None):
    """ write a series of geometries to an .xyz trajectory file
    """
    ngeos = len(xyzs_lst)
    comments = ('',) * ngeos if comments is None else comments
    assert len(comments) == ngeos
    xyz_strs = [write_xyz(syms, xyzs, comment)
                for xyzs, comment in zip(xyzs_lst, comments)]
    return '\n'.join(xyz_strs)
