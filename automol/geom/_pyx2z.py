""" pyx2z interface
"""

import pyx2z
import autoread as ar
import autoparse.pattern as app
import automol.zmat.base


def to_oriented_geometry(geo):
    """ Generate an oriented x2z molecule object from a geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: x2z molecule object
    """

    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)
    orient_mg = pyx2z.MolecOrient(_mg)

    return orient_mg


def from_geometry(geo, ts_bnds=()):
    """ Generate an x2z molecule object from a geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param ts_bnds: keys for the breaking/forming bonds in a TS
        :type ts_bnds: tuple(frozenset(int))
        :rtype: x2z molecule object
    """

    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)

    ts_bnds = list(map(list, ts_bnds))
    x2m = pyx2z.MolecStruct(_mg, ts_bnds)

    return x2m


def to_zmatrix(x2m):
    """ Generate a Z-Natrix from an x2z molecule object.

        :param x2m: molecule object
        :type x2m: x2z molecule object
        :rtype: automol Z-Matrix data structure
    """

    zma_str = pyx2z.zmatrix_string(x2m)
    syms, key_mat, name_mat, val_mat = ar.zmat.read(
        zma_str,
        mat_entry_sep_ptt=',',
        mat_entry_start_ptt=',',
        setv_sep_ptt=app.padded(app.one_of_these(['', app.NEWLINE])))

    zma = automol.zmat.base.from_data(
        syms, key_mat, val_mat, name_mat,
        one_indexed=True, angstrom=False, degree=True)

    return zma


def zmatrix_torsion_coordinate_names(x2m):
    """ Build a list of Z-Matrix torsion coordinate names
        from an x2z molecule object.

        :param x2m: molecule object
        :type x2m: x2z molecule object
        :rtype: tuple(str)
    """
    return pyx2z.rotational_bond_coordinates(x2m)


def zmatrix_atom_ordering(x2m):
    """ Build dictionary that acts as a mapping from the order of atoms in the
        x2z molecule object to the order of atoms in the reslutant Z-Matrix.

        :param x2m: molecule object
        :type x2m: x2z molecule object
        :rtype: dict[int: int]
    """
    return {geo_key: zma_key
            for zma_key, geo_key in enumerate(x2m.atom_ordering())}
