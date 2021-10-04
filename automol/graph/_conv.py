""" graph conversions
"""

import autoparse.pattern as app
import autoparse.find as apf
from automol.util import dict_
import automol.graph.embed
import automol.geom.base
import automol.inchi.base
from automol.extern import molfile
from automol.extern import rdkit_
from automol.graph.base import atom_keys
from automol.graph.base import bond_keys
from automol.graph.base import atom_symbols
from automol.graph.base import bond_orders
from automol.graph.base import atom_bond_valences
from automol.graph.base import atom_unsaturated_valences
from automol.graph.base import has_stereo
from automol.graph.base import explicit
from automol.graph.base import without_dummy_atoms
from automol.graph.base import backbone_isomorphic
from automol.graph.base import dominant_resonance


# # conversions
def geometry(gra):
    """ Convert a molecular graph to a molecular geometry.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :rtype: automol molecular geometry data structure
    """

    symbs = atom_symbols(gra)
    if len(symbs) != 1:
        gra = explicit(gra)
        geo = automol.graph.embed.geometry(gra)
    else:
        symb = list(symbs.values())[0]
        # symb = list(symbs.keys())[0]
        geo = ((symb, (0.00, 0.00, 0.00)),)

    return geo


def inchi(gra, stereo=False):
    """ Generate an InChI string from a molecular graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """

    ich = automol.inchi.base.hardcoded_object_to_inchi_by_key(
        'graph', gra, comp=_compare)

    if ich is None:
        if not stereo or not has_stereo(gra):
            ich, _ = inchi_with_sort_from_geometry(gra)
            ich = automol.inchi.base.standard_form(ich, stereo=stereo)
        else:
            gra = explicit(gra)
            geo, geo_idx_dct = automol.graph.embed.fake_stereo_geometry(gra)
            ich, _ = inchi_with_sort_from_geometry(
                gra, geo=geo, geo_idx_dct=geo_idx_dct)

    return ich


def stereo_inchi(gra):
    """ Generate an InChI string from a molecular graph, including stereo
        information.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :rtype: str
    """
    return inchi(gra, stereo=True)


def inchi_with_sort_from_geometry(gra, geo=None, geo_idx_dct=None):
    """ Generate an InChI string from a molecular graph.
        If coordinates are passed in, they are used to determine stereo.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct:
        :type geo_idx_dct: dict[:]
        :rtype: (str, tuple(int))
    """
    mlf, key_map_inv = molfile_with_atom_mapping(gra, geo=geo,
                                                 geo_idx_dct=geo_idx_dct)
    rdm = rdkit_.from_molfile(mlf)
    ich, aux_info = rdkit_.to_inchi(rdm, with_aux_info=True)
    nums = _parse_sort_order_from_aux_info(aux_info)
    nums = tuple(map(key_map_inv.__getitem__, nums))

    return ich, nums


def molfile_with_atom_mapping(gra, geo=None, geo_idx_dct=None):
    """ Generate an MOLFile from a molecular graph.
        If coordinates are passed in, they are used to determine stereo.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct:
        :type geo_idx_dct: dict[:]
        :returns: the MOLFile string, followed by a mapping from MOLFile atoms
            to atoms in the graph
        :rtype: (str, dict)
    """
    gra = without_dummy_atoms(gra)
    gra = dominant_resonance(gra)
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = list(bond_keys(gra))
    atm_syms = dict_.values_by_key(atom_symbols(gra), atm_keys)
    atm_bnd_vlcs = dict_.values_by_key(
        atom_bond_valences(gra), atm_keys)
    atm_rad_vlcs = dict_.values_by_key(
        atom_unsaturated_valences(gra), atm_keys)
    bnd_ords = dict_.values_by_key(bond_orders(gra), bnd_keys)

    if geo is not None:
        assert geo_idx_dct is not None
        atm_xyzs = automol.geom.base.coordinates(geo)
        atm_xyzs = [atm_xyzs[geo_idx_dct[atm_key]] if atm_key in geo_idx_dct
                    else (0., 0., 0.) for atm_key in atm_keys]
    else:
        atm_xyzs = None

    mlf, key_map_inv = molfile.from_data(
        atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs, atm_rad_vlcs, bnd_ords,
        atm_xyzs=atm_xyzs)
    return mlf, key_map_inv


def rdkit_molecule(gra):
    """ Convert a molecular graph to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(gra))

    :param gra: the graph
    :returns: the RDKit molecule
    """
    return rdkit_.from_inchi(inchi(gra))


# # helpers
def _parse_sort_order_from_aux_info(aux_info):
    ptt = app.escape('/N:') + app.capturing(
        app.series(app.UNSIGNED_INTEGER, ','))
    num_str = apf.first_capture(ptt, aux_info)
    nums = tuple(map(int, num_str.split(',')))
    return nums


def _compare(gra1, gra2):
    """ Compare the backbone structure of two moleculare graphs.

        :param gra1: molecular graph 1
        :type gra1: automol graph data structure
        :param gra2: molecular graph 2
        :type gra2: automol graph data structure
        :rtype: bool
    """

    gra1 = without_dummy_atoms(gra1)
    gra2 = without_dummy_atoms(gra2)

    return backbone_isomorphic(gra1, gra2)
