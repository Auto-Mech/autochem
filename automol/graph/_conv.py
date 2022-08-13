""" graph conversions
"""
import functools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from automol import error
import automol.graph.embed
import automol.geom.base
import automol.inchi.base
from automol.util import dict_
from automol.extern import molfile
from automol.extern import rdkit_
from automol.extern import pybel_
from automol.graph.base import atom_keys
from automol.graph.base import bond_keys
from automol.graph.base import atom_symbols
from automol.graph.base import bond_orders
from automol.graph.base import atom_bond_valences
from automol.graph.base import atom_unsaturations
from automol.graph.base import explicit
from automol.graph.base import without_dummy_atoms
from automol.graph.base import backbone_isomorphic
from automol.graph.base import kekule
from automol.graph.base import set_stereo_from_geometry
from automol.graph.base import smiles
from automol.graph.base import amchi
from automol.graph.base import inchi_is_bad
from automol.graph.base import implicit
from automol.graph.base import to_local_stereo
from automol.graph.base import subgraph
from automol.graph.base import relabel
from automol.graph.base import bond_stereo_keys
from automol.graph.base import explicit_hydrogen_keys
from automol.graph.base import bond_stereo_parities
from automol.graph.base import connected_components
from automol.graph.base import has_stereo
from automol.graph.base import string


# # conversions
def geometry(gra, check=True):
    """ Convert a molecular graph to a molecular geometry.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry data structure
    """
    gra = automol.graph.explicit(gra)
    gras = connected_components(gra)

    geos = [_connected_geometry(g, check=check) for g in gras]
    geos = [automol.geom.translate(g, [50.*i, 0., 0.])
            for i, g in enumerate(geos)]
    geo = functools.reduce(automol.geom.join, geos)

    # Work out how the geometry needs to be re-ordered to match the input graph
    idx_dct = {}
    idxs_lst = [sorted(atom_keys(g)) for g in gras]
    offset = 0
    for idxs in idxs_lst:
        idx_dct.update({i+offset: k for i, k in enumerate(idxs)})
        offset += len(idxs)

    geo = automol.geom.reorder(geo, idx_dct)
    return geo


def _connected_geometry(gra, check=True):
    """ Generate a geometry for a connected molecular graph.

        :param gra: connected molecular graph
        :type gra: automol graph data structure
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry
    """
    ste_keys = automol.graph.stereo_keys(gra)

    smi = smiles(gra, res_stereo=False)
    has_ste = has_stereo(gra)

    def _gen1():
        nonlocal smi

        rdm = rdkit_.from_smiles(smi)
        geo, = rdkit_.to_conformers(rdm, nconfs=1)
        return geo

    def _gen2():
        nonlocal smi

        pbm = pybel_.from_smiles(smi)
        geo = pybel_.to_geometry(pbm)
        return geo

    def _gen3():
        nonlocal gra

        if has_ste:
            raise ValueError

        gra_ = automol.graph.explicit(gra)
        geo = automol.graph.embed.geometry(gra_)
        return geo

    success = False
    for gen_ in (_gen1, _gen1, _gen1, _gen2, _gen3):
        try:
            geo = gen_()
        except (RuntimeError, TypeError, ValueError):
            continue

        if check:
            # First, check connectivity.
            gra_ = automol.geom.graph(geo)
            geo = automol.graph.linear_vinyl_corrected_geometry(gra_, geo)

            idx_dct = automol.graph.isomorphism(gra_, gra, stereo=False)

            if idx_dct is None:
                continue

            # Reorder the geometry to match the input graph connectivity.
            geo = automol.geom.reorder(geo, idx_dct)

            # If connectivity matches and there is no stereo, we are done.
            if not has_ste:
                success = True
                break

            # Otherwise, there is stereo.
            # First, try an isomorphism to see if the parities already match.
            gra_ = automol.geom.graph(geo)
            par_dct_ = automol.graph.stereo_parities(gra_)
            par_dct_ = {k: (p if k in ste_keys else None)
                        for k, p in par_dct_.items()}
            gra_ = automol.graph.set_stereo_parities(gra_, par_dct_)
            idx_dct = automol.graph.isomorphism(gra_, gra, stereo=True)

            # If connectivity and stereo match, we are done.
            if idx_dct:
                # Reorder the geometry to match the input graph
                geo = automol.geom.reorder(geo, idx_dct)
                success = True
                break

            # If the stereo doesn't match, try a stereo correction.
            geo = automol.graph.stereo_corrected_geometry(
                gra, geo, local_stereo=False)

            # Now, re-try the isomorphism
            gra_ = automol.geom.graph(geo)
            par_dct_ = automol.graph.stereo_parities(gra_)
            par_dct_ = {k: (p if k in ste_keys else None)
                        for k, p in par_dct_.items()}
            gra_ = automol.graph.set_stereo_parities(gra_, par_dct_)
            idx_dct = automol.graph.isomorphism(gra_, gra, stereo=True)

            # If this fails, this geometry won't work. Continue
            if not idx_dct:
                continue

            # The stereo matches after correction.
            # Reorder the geometry to match the input graph, and we are done.
            geo = automol.geom.reorder(geo, idx_dct)
            success = True
            break

    if not success:
        raise error.FailedGeometryGenerationError(
            f'Failed graph:\n{string(gra)}')

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

    # Assuming the MolFile InChI works, the above code is all we need. What
    # follows is to correct cases where it fails.
    # This only appears to work sometimes, so when it doesn't, we fall back on
    # the original inchi output.
    if geo is not None:
        gra = set_stereo_from_geometry(gra, geo, geo_idx_dct=geo_idx_dct)
        gra = implicit(gra)
        sub_ichs = automol.inchi.split(ich)

        failed = False

        new_sub_ichs = []
        for sub_ich, nums in zip(sub_ichs, nums_lst):
            sub_gra = subgraph(gra, nums, stereo=True)
            sub_ich = _connected_inchi_with_graph_stereo(
                sub_ich, sub_gra, nums)
            if sub_ich is None:
                failed = True
                break

            new_sub_ichs.append(sub_ich)

        # If it worked, replace the InChI with our forced-stereo InChI.
        if not failed:
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
        # Convert to implicit graph and relabel based on InChI sort
        atm_key_dct = dict(map(reversed, enumerate(nums)))
        gra = relabel(gra, atm_key_dct)
        gra = explicit(gra)
        exp_h_keys = explicit_hydrogen_keys(gra)
        exp_h_key_dct = {k: -k for k in exp_h_keys}
        gra = relabel(gra, exp_h_key_dct)

        gra = to_local_stereo(gra)

        # Translate internal stereo parities into InChI stereo parities
        # and generate the appropriate b-layer string for the InChI
        ste_dct = bond_stereo_parities(gra)
        ste_keys = tuple(sorted(tuple(reversed(sorted(k)))
                                for k in bond_stereo_keys(gra)))
        blyr_strs = []
        for atm1_key, atm2_key in ste_keys:
            par = ste_dct[frozenset({atm1_key, atm2_key})]

            blyr_strs.append(
                f"{atm1_key+1}-{atm2_key+1}{'+' if par else '-'}")

        # After forming the b-layer string, generate the new InChI
        blyr_str = ','.join(blyr_strs)
        ste_dct = {'b': blyr_str}
        # print(ste_dct)
        ich = automol.inchi.standard_form(ich, ste_dct=ste_dct)
        # print('out:', ich)

    return ich


def chi(gra, stereo=True):
    """ Generate a ChI string from a molecular graph.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :returns: ChI string
        :rtype: str
    """
    ret = inchi(gra, stereo=stereo)
    if inchi_is_bad(gra, ret):
        ret = amchi(gra, stereo=stereo)

    return ret


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
    gra = kekule(gra)
    atm_keys = sorted(atom_keys(gra))
    bnd_keys = list(bond_keys(gra))
    atm_syms = dict_.values_by_key(atom_symbols(gra), atm_keys)
    atm_bnd_vlcs = dict_.values_by_key(
        atom_bond_valences(gra), atm_keys)
    atm_rad_vlcs = dict_.values_by_key(
        atom_unsaturations(gra), atm_keys)
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
