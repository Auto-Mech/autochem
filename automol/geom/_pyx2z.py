""" pyx2z interface
"""
import pyx2z
from automol.constructors.zmatrix import from_data as _zmatrix_from_data
import autoread as ar
import autoparse.pattern as app


def from_geometry(geo):
    """ x2z molecule object from a geometry
    """
    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)
    _ps = pyx2z.PrimStruct(_mg)
    x2m = pyx2z.MolecStruct(_ps)
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

    zma = _zmatrix_from_data(
        syms, key_mat, name_mat, val_dct,
        one_indexed=True, angstrom=False, degree=True, complete=True)
    return zma


def zmatrix_torsion_coordinate_names(x2m):
    """ z-matrix torsion coordinate name from an x2z molecule object
    """
    names = pyx2z.rotational_bond_coordinates(x2m)
    return names
