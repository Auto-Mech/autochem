""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.inchi.base
import automol.graph.base
from automol.extern import rdkit_
from automol.smiles.base import split
from automol.smiles.base import without_resonance_stereo
from automol.smiles.base import parse_connected_molecule_properties


# # conversions
def amchi(smi):
    """ Generate an AMChI string from a connected SMILES string.

        :param smi: SMILES string
        :type smi: str
        :returns: AMChI string
        :rtype: str
    """
    gra = graph(smi, stereo=True, local_stereo=False)
    ach = automol.graph.base.amchi(gra)
    return ach


def inchi(smi):
    """ Convert a SMILES string into an InChI string.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """
    smi = without_resonance_stereo(smi)

    ich = automol.inchi.base.hardcoded_object_to_inchi_by_key(
        'smiles', smi, comp=_compare)

    if ich is None:
        rdm = rdkit_.from_smiles(smi)
        ich = rdkit_.to_inchi(rdm)
    return ich


def chi(smi):
    """ Convert a SMILES string to an AMChI or InChI string.

        Currently only uses AMChI for resonance bond stereo.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """
    gra = graph(smi, stereo=True, local_stereo=False)
    ret = inchi(smi)
    if automol.graph.base.inchi_is_bad(gra, ret):
        ret = automol.graph.base.amchi(gra)

    return ret


def graph(smi, stereo=True, local_stereo=False):
    """ Generate a molecular graph from a SMILES string.

        :param smi: SMILES string
        :type smi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :param local_stereo: assign local stereo parities?
        :type local_stereo: bool
        :rtype: automol molecular graph
    """
    smis = split(smi)
    gras = [_connected_graph(s, stereo=stereo, local_stereo=local_stereo)
            for s in smis]
    gra = automol.graph.base.union_from_sequence(gras, shift_keys=True)
    return gra


def _connected_graph(smi, stereo=True, local_stereo=False):
    """ Generate a connected molecular graph from a connected SMILES string.

        :param smi: SMILES string
        :type smi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :param local_stereo: assign local stereo parities?
        :type local_stereo: bool
        :rtype: automol molecular graph
    """
    symb_dct, bnd_ord_dct, atm_par_dct, bnd_par_dct = (
        parse_connected_molecule_properties(smi))
    bnd_keys = bnd_ord_dct.keys()

    if not stereo:
        atm_par_dct = None
        bnd_par_dct = None

    gra = automol.graph.base.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_ste_par_dct=atm_par_dct,
        bnd_ste_par_dct=bnd_par_dct,
    )

    if automol.graph.base.has_stereo(gra):
        # The parser marks all bonds with directional bonds on either side as
        # having stereo, because it has no way to distinguish between them.  In
        # lieu of a more rigorous check, remove stereo from all non-sp2 bonds.
        # If this is an issue, we could create a more rigorous check to see if
        # a bond is stereogenic.
        ste_bnd_keys = automol.graph.base.bond_stereo_keys(gra)
        sp2_bnd_keys = automol.graph.base.sp2_bond_keys(gra)
        bnd_keys = ste_bnd_keys - sp2_bnd_keys
        gra = automol.graph.base.remove_bond_stereo_parities(gra, bnd_keys)

        if not local_stereo:
            # Convert from local to canonical stereo
            gra = automol.graph.base.from_local_stereo(gra)

    return gra


# helpers
def _compare(smi1, smi2):
    """ Check if two SMILES strings are similar.

        :param smi1: SMILES string 1
        :type smi1: str
        :param smi2: SMILES string 2
        :type smi2: str
        :rtype: bool
    """
    return _canonicalize(smi1) == _canonicalize(smi2)


def _canonicalize(smi):
    """ Convert a SMILES string into its canonical form.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """
    return rdkit_.to_smiles(rdkit_.from_smiles(smi))


if __name__ == '__main__':
    SMIS = [
        '[CH2]CC[C@H]1O[C@H]1C.O',
        '[CH2]CC[C@@H]1O[C@@H]1C.[OH]',
    ]
    for SMI in SMIS:
        ACH = amchi(SMI)
        print(ACH)
