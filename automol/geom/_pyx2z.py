""" pyx2z interface
"""
import importlib.util

from automol.zmat import base as zmat_base

pyx2z_found = importlib.util.find_spec("pyx2z")
if pyx2z_found is not None:
    import pyx2z
else:
    pyx2z = None

autoparse_found = importlib.util.find_spec("autoparse")
if autoparse_found is not None:
    import autoparse.pattern as app
else:
    app = None

autoread_found = importlib.util.find_spec("autoread")
if autoread_found is not None:
    import autoread as ar
else:
    ar = None


def from_geometry(geo, ts_bnds=()):
    """Generate an x2z molecule object from a geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ts_bnds: keys for the breaking/forming bonds in a TS
    :type ts_bnds: tuple(frozenset(int))
    :rtype: x2z molecule object
    """
    if pyx2z is None:
        raise NotImplementedError("Not implemented without x2z!")

    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)

    ts_bnds = list(map(list, ts_bnds))
    x2m = pyx2z.MolecStruct(_mg, ts_bnds)

    return x2m


def to_zmatrix(x2m):
    """Generate a Z-Natrix from an x2z molecule object.

    :param x2m: molecule object
    :type x2m: x2z molecule object
    :rtype: automol Z-Matrix data structure
    """
    if pyx2z is None:
        raise NotImplementedError("Not implemented without x2z!")

    if ar is None:
        raise NotImplementedError("Not implemented without autoread!")

    if app is None:
        raise NotImplementedError("Not implemented without autoparse!")

    zma_str = pyx2z.zmatrix_string(x2m)
    syms, key_mat, name_mat, val_mat = ar.zmat.read(
        zma_str,
        mat_entry_sep_ptt=",",
        mat_entry_start_ptt=",",
        setv_sep_ptt=app.padded(app.one_of_these(["", app.NEWLINE])),
    )

    zma = zmat_base.from_data(
        syms, key_mat, val_mat, name_mat, one_indexed=True, angstrom=False, degree=True
    )

    return zma


def zmatrix_torsion_coordinate_names(x2m):
    """Build a list of Z-Matrix torsion coordinate names
    from an x2z molecule object.

    :param x2m: molecule object
    :type x2m: x2z molecule object
    :rtype: tuple(str)
    """
    if pyx2z is None:
        raise NotImplementedError("Not implemented without x2z!")

    return pyx2z.rotational_bond_coordinates(x2m)
