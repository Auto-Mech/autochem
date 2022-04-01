""" graph conversions
"""

import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
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
from automol.graph.base import explicit
from automol.graph.base import without_dummy_atoms
from automol.graph.base import backbone_isomorphic
from automol.graph.base import dominant_resonance
from automol.graph.base import set_stereo_from_geometry
from automol.graph.base import smiles


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


def inchi(gra, stereo=True):
    """ Generate an InChI string from a molecular graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """
    smi = smiles(gra, stereo=stereo, res_stereo=False)

    ich = automol.inchi.base.hardcoded_object_to_inchi_by_key(
        'smiles', smi, comp=_compare_smiles)

    if ich is None:
        rdm = rdkit_.from_smiles(smi)
        ich = rdkit_.to_inchi(rdm)
    return ich


def inchi_with_sort_from_geometry(gra, geo=None, geo_idx_dct=None):
    """ Generate an InChI string from a molecular graph.
        If coordinates are passed in, they are used to determine stereo.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct:
        :type geo_idx_dct: dict[:]
        :returns: the inchi string, along with the InChI sort order of the
            atoms
        :rtype: (str, tuple(int))
    """
    if geo is not None:
        natms = automol.geom.base.count(geo)
        geo_idx_dct = (dict(enumerate(range(natms)))
                       if geo_idx_dct is None else geo_idx_dct)
        gra = set_stereo_from_geometry(gra, geo, geo_idx_dct=geo_idx_dct)

    mlf, key_map_inv = molfile_with_atom_mapping(gra, geo=geo,
                                                 geo_idx_dct=geo_idx_dct)
    rdm = rdkit_.from_molfile(mlf)
    ich, aux_info = rdkit_.to_inchi(rdm, with_aux_info=True)

    nums_lst = _parse_sort_order_from_aux_info(aux_info)
    nums_lst = tuple(tuple(map(key_map_inv.__getitem__, nums))
                     for nums in nums_lst)

    return ich, nums_lst


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
        app.series(app.UNSIGNED_INTEGER, app.one_of_these(',;')))
    num_strs = apf.first_capture(ptt, aux_info).split(';')
    nums_lst = ap_cast(tuple(s.split(',') for s in num_strs))
    return nums_lst


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


def _compare_smiles(smi1, smi2):
    """ Check if two SMILES strings are similar.

        :param smi1: SMILES string 1
        :type smi1: str
        :param smi2: SMILES string 2
        :type smi2: str
        :rtype: bool
    """
    return _canonicalize_smiles(smi1) == _canonicalize_smiles(smi2)


def _canonicalize_smiles(smi):
    """ Convert a SMILES string into its canonical form.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """
    return rdkit_.to_smiles(rdkit_.from_smiles(smi))
