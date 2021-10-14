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
from automol.graph.base import has_stereo
from automol.graph.base import explicit
from automol.graph.base import implicit
from automol.graph.base import relabel
from automol.graph.base import subgraph
from automol.graph.base import without_dummy_atoms
from automol.graph.base import backbone_isomorphic
from automol.graph.base import dominant_resonance
from automol.graph.base import bond_stereo_keys
from automol.graph.base import bond_stereo_parities
from automol.graph.base import bond_stereo_sorted_neighbor_atom_keys
from automol.graph.base import set_stereo_from_geometry


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
            ich, nums = inchi_with_sort_from_geometry(
                gra, geo=geo, geo_idx_dct=geo_idx_dct)

            # First, do a check to see if the InChI is missing stereo relative
            # to the graph.
            # >>> check here

            # Convert to an implicit graph and relabel based on the InChI sort
            gra = implicit(gra)
            atm_key_dct = dict(map(reversed, enumerate(nums)))
            print(atm_key_dct)
            print(nums)
            gra = relabel(gra, atm_key_dct)
            print(automol.graph.string(gra))

            ste_dct = bond_stereo_parities(gra)
            ste_keys = sorted(map(sorted, bond_stereo_keys(gra)))
            ste_vals = dict_.values_by_key(ste_dct, map(frozenset, ste_keys))
            print(ste_keys)
            print(ste_vals)

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
        :returns: the inchi string, along with the InChI sort order of the
            atoms
        :rtype: (str, tuple(int))
    """
    mlf, key_map_inv = molfile_with_atom_mapping(gra, geo=geo,
                                                 geo_idx_dct=geo_idx_dct)
    rdm = rdkit_.from_molfile(mlf)
    ich, aux_info = rdkit_.to_inchi(rdm, with_aux_info=True)

    nums_lst = _parse_sort_order_from_aux_info(aux_info)
    nums_lst = tuple(tuple(map(key_map_inv.__getitem__, nums))
                     for nums in nums_lst)

    # Assuming the MolFile InChI works, the above code is all we need. What
    # follows is to correct cases where it fails.
    # Limitation: This code could fail in cases where the InChI string is
    # missing two double bonds that are isomorphically equivalent.
    if geo is not None:
        gra = set_stereo_from_geometry(gra, geo, geo_idx_dct=geo_idx_dct)
        gra = implicit(gra)
        sub_ichs = automol.inchi.split(ich)

        new_sub_ichs = []
        for sub_ich, nums in zip(sub_ichs, nums_lst):
            sub_gra = subgraph(gra, nums, stereo=True)
            sub_ich = _connected_inchi_with_graph_stereo(
                sub_ich, sub_gra, nums)
            new_sub_ichs.append(sub_ich)

        ich = automol.inchi.join(new_sub_ichs)
        ich = automol.inchi.standard_form(ich)

    return ich, nums_lst


def _connected_inchi_with_graph_stereo(ich, gra, nums):
    """ For a connected inchi/graph, check if the inchi is missing stereo; If so,
    add stereo based on the graph.

    Currently only checks for missing bond stereo, since this is all we have
    seen so far, but could be generalized.

    :param ich: the inchi string
    :param gra: the graph
    :param nums: graph indices to backbone atoms in canonical inchi order
    :type nums: tuple[int]
    """
    # First, do a check to see if the InChI is missing bond stereo
    # relative to the graph.
    ich_ste_keys = automol.inchi.stereo_bonds(ich)
    our_ste_keys = bond_stereo_keys(gra)

    miss_ich_ste_keys = automol.inchi.unassigned_stereo_bonds(ich)

    if len(ich_ste_keys) > len(our_ste_keys):
        raise Exception("Our code is missing stereo bonds")

    if len(ich_ste_keys) < len(our_ste_keys) or miss_ich_ste_keys:
        print(automol.graph.string(gra))
        # Convert to implicit graph and relabel based on InChI sort
        atm_key_dct = dict(map(reversed, enumerate(nums)))
        gra = relabel(gra, atm_key_dct)

        # Translate internal stereo parities into InChI stereo parities
        # and generate the appropriate b-layer string for the InChI
        ste_dct = bond_stereo_parities(gra)
        ste_keys = tuple(sorted(tuple(reversed(sorted(k)))
                                for k in bond_stereo_keys(gra)))
        blyr_strs = []
        for atm1_key, atm2_key in ste_keys:
            our_par = ste_dct[frozenset({atm1_key, atm2_key})]
            our_srt1, our_srt2 = bond_stereo_sorted_neighbor_atom_keys(
                gra, atm1_key, atm2_key)
            ich_srt1 = tuple(reversed(our_srt1))
            ich_srt2 = tuple(reversed(our_srt2))
            if not ((our_srt1 != ich_srt1) ^ (our_srt2 != ich_srt2)):
                ich_par = our_par
            else:
                ich_par = not our_par

            blyr_strs.append(
                f"{atm1_key+1}-{atm2_key+1}{'-' if ich_par else '+'}")

        # After forming the b-layer string, generate the new InChI
        blyr_str = ','.join(blyr_strs)
        ste_dct = {'b': blyr_str}
        ich = automol.inchi.standard_form(ich, ste_dct=ste_dct)

    return ich


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
