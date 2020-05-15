""" geometry conversions
"""
import itertools
import numpy
from automol import create
from automol.convert import _pyx2z
from automol.convert import _util
import automol.graph
import automol.geom
import automol.zmatrix
import automol.convert.graph
import automol.convert.inchi


# geometry => z-matrix
def zmatrix(geo, ts_bnds=()):
    """ geometry => z-matrix
    """
    syms = automol.geom.symbols(geo)
    if len(syms) == 1:
        key_mat = [[None, None, None]]
        name_mat = [[None, None, None]]
        val_dct = {}
        zma = create.zmatrix.from_data(syms, key_mat, name_mat, val_dct)
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        zma = _pyx2z.to_zmatrix(x2m)
    zma = automol.zmatrix.standard_form(zma)
    return zma


def zmatrix_torsion_coordinate_names(geo, ts_bnds=()):
    """ z-matrix torsional coordinate names
    """
    syms = automol.geom.symbols(geo)
    if len(syms) == 1:
        names = ()
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        names = _pyx2z.zmatrix_torsion_coordinate_names(x2m)

        zma = _pyx2z.to_zmatrix(x2m)
        name_dct = automol.zmatrix.standard_names(zma)
        names = tuple(map(name_dct.__getitem__, names))
    return names


def zmatrix_atom_ordering(geo, ts_bnds=()):
    """ z-matrix atom ordering
    """
    syms = automol.geom.symbols(geo)
    if len(syms) == 1:
        idxs = (0,)
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        idxs = _pyx2z.zmatrix_atom_ordering(x2m)
    return idxs


# geometry => graph
def graph(geo, remove_stereo=False):
    """ geometry => graph
    """
    gra = _connectivity_graph(geo)
    if not remove_stereo:
        gra = automol.graph.set_stereo_from_geometry(gra, geo)
    return gra


def weakly_connected_graph(geo, remove_stereo=False):
    """ geometry => graph
    """
    gra = _connectivity_graph(
        geo, rqq_bond_max=5.0, rqh_bond_max=5.0, rhh_bond_max=2.3)
    _ = remove_stereo
    # if not remove_stereo:
    #    xyzs = automol.geom.coordinates(geo)
    #    atm_xyz_dct = dict(enumerate(xyzs))
    #    gra = automol.graph.set_stereo_from_atom_coordinates(gra, atm_xyz_dct)
    return gra


def _connectivity_graph(geo,
                        rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ geometry => connectivity graph (no stereo)
    """
    syms = automol.geom.symbols(geo)
    xyzs = automol.geom.coordinates(geo)

    def _are_bonded(idx_pair):
        xyz1, xyz2 = map(xyzs.__getitem__, idx_pair)
        sym1, sym2 = map(syms.__getitem__, idx_pair)
        dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))
        return (False if 'X' in (sym1, sym2) else
                (dist < rqh_bond_max) if 'H' in (sym1, sym2) else
                (dist < rhh_bond_max) if (sym1 == 'H' and sym2 == 'H') else
                (dist < rqq_bond_max))

    idxs = range(len(xyzs))
    atm_sym_dct = dict(enumerate(syms))
    bnd_keys = tuple(
        map(frozenset, filter(_are_bonded, itertools.combinations(idxs, r=2))))
    gra = create.graph.from_data(atom_symbols=atm_sym_dct, bond_keys=bnd_keys)
    return gra


# geometry => inchi
def inchi(geo, remove_stereo=False):
    """ geometry => InChI
    """
    ich, _ = inchi_with_sort(geo, remove_stereo=remove_stereo)
    return ich


def inchi_with_sort(geo, remove_stereo=False):
    """ geometry => InChI
    """
    ich = automol.convert.inchi.object_to_hardcoded_inchi_by_key(
        'geom', geo, comp=_compare)
    nums = None

    if ich is None:
        gra = _connectivity_graph(geo)
        if remove_stereo:
            geo = None
            geo_idx_dct = None
        else:
            geo_idx_dct = dict(enumerate(range(automol.geom.count(geo))))
        ich, nums = automol.convert.graph.inchi_with_sort_from_geometry(
            gra=gra, geo=geo, geo_idx_dct=geo_idx_dct)

    return ich, nums


def _compare(geo1, geo2):
    gra1 = automol.graph.without_dummy_atoms(_connectivity_graph(geo1))
    gra2 = automol.graph.without_dummy_atoms(_connectivity_graph(geo2))
    return automol.graph.backbone_isomorphic(gra1, gra2)


# geometry => formula
def formula(geo):
    """ geometry => formula
    """
    syms = automol.geom.symbols(geo)
    fml = _util.formula(syms)
    return fml
