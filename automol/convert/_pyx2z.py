""" pyx2z interface
"""
import numpy
import pyx2z
import autoread as ar
import autoparse.pattern as app
import automol.create.zmat


def to_oriented_geometry(geo):
    """ Obtain an oriented x2z molecule object from a geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :type: x2z molecule object
    """

    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)
    orient_mg = pyx2z.MolecOrient(_mg)
    return orient_mg


def from_geometry(geo, ts_bnds=()):
    """ x2z molecule object from a geometry
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
    """ z-matrix from an x2z molecule object
    """
    zma_str = pyx2z.zmatrix_string(x2m)

    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(
        zma_str,
        mat_entry_sep_ptt=',',
        mat_entry_start_ptt=',',
        setv_sep_ptt=app.padded(app.one_of_these(['', app.NEWLINE])))

    val_dct[None] = None
    val_mat = numpy.vectorize(val_dct.__getitem__)(name_mat)
    zma = automol.create.zmat.from_data(
        syms, key_mat, val_mat, name_mat,
        one_indexed=True, angstrom=False, degree=True)
    return zma


def to_zmatrix_old(x2m):
    """ z-matrix from an x2z molecule object
    """
    zma_str = pyx2z.zmatrix_string(x2m)

    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(
        zma_str,
        mat_entry_sep_ptt=',',
        mat_entry_start_ptt=',',
        setv_sep_ptt=app.padded(app.one_of_these(['', app.NEWLINE])))

    zma = automol.create.zmatrix.from_data(
        syms, key_mat, name_mat, val_dct,
        one_indexed=True, angstrom=False, degree=True)

    return zma


def zmatrix_torsion_coordinate_names(x2m):
    """ z-matrix torsion coordinate name from an x2z molecule object
    """
    names = pyx2z.rotational_bond_coordinates(x2m)
    return names


def zmatrix_atom_ordering(x2m):
    """ z-matrix atom ordering from an x2z molecule object
    """
    idx_dct = {geo_key: zma_key
               for zma_key, geo_key in enumerate(x2m.atom_ordering())}
    return idx_dct
