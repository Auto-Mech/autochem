""" working with geometries

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import numpy
from phydat import phycon
from automol import util
import automol.geom.base
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bonds_neighbor_atom_keys
from automol.graph.base._core import bonds_neighbor_bond_keys
from automol.graph.base._algo import branch
from automol.graph.base._resonance import sp2_bond_keys


def linear_vinyl_corrected_geometry(gra, geo, geo_idx_dct=None,
                                    tol=2.*phycon.DEG2RAD):
    """ correct a geometry for linear vinyl groups

        :param gra: molecular graph with stereo parities
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :param tol: tolerance of bond angle(s) for determing linearity
        :type tol: float
    """
    atm_keys = atom_keys(gra)
    bnakeys_dct = bonds_neighbor_atom_keys(gra)
    bnbkeys_dct = bonds_neighbor_bond_keys(gra)

    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(sorted(atm_keys))})

    bnd_keys = sp2_bond_keys(gra)

    for bnd1_key in bnd_keys:
        for bnd2_key in bnbkeys_dct[bnd1_key]:
            atm2_key, = bnd1_key & bnd2_key
            atm1_key, = bnd1_key - {atm2_key}
            atm3_key, = bnd2_key - {atm2_key}

            atm1_idx = geo_idx_dct[atm1_key]
            atm2_idx = geo_idx_dct[atm2_key]
            atm3_idx = geo_idx_dct[atm3_key]

            ang = automol.geom.base.central_angle(
                geo, atm1_idx, atm2_idx, atm3_idx)

            if numpy.abs(ang - numpy.pi) < tol:
                atm0_key = next(iter(bnakeys_dct[bnd1_key] - {atm3_key}), None)
                atm0_key = atm3_key if atm0_key is None else atm0_key
                atm0_idx = geo_idx_dct[atm0_key]

                xyzs = automol.geom.base.coordinates(geo)

                atm0_xyz = xyzs[atm0_idx]
                atm1_xyz = xyzs[atm1_idx]
                atm2_xyz = xyzs[atm2_idx]

                rot_axis = util.vec.unit_perpendicular(atm0_xyz, atm1_xyz,
                                                       orig_xyz=atm2_xyz)

                rot_atm_keys = atom_keys(
                    branch(gra, atm2_key, {atm2_key, atm3_key}))

                rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

                geo = automol.geom.rotate(
                    geo, rot_axis, numpy.pi/3,
                    orig_xyz=atm2_xyz, idxs=rot_idxs)

    return geo
