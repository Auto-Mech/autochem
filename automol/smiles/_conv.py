""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.inchi.base
import automol.graph.base
from automol.extern import rdkit_
from automol.smiles.base import split
from automol.smiles.base import without_resonance_stereo
from automol.smiles.base import without_stereo
from automol.smiles.base import parse_connected_molecule_properties


# # conversions
def amchi(smi, stereo=True):
    """Generate an AMChI string from a connected SMILES string.

    :param smi: SMILES string
    :type smi: str
    :param stereo: Keep stereo information form the SMILES string?
    :type stereo: bool
    :returns: AMChI string
    :rtype: str
    """
    gra = graph(smi, stereo=stereo, local_stereo=False)
    ach = automol.graph.base.amchi(gra)
    return ach


def inchi(smi, stereo=True):
    """Convert a SMILES string into an InChI string.

    :param smi: SMILES string
    :type smi: str
    :param stereo: Keep stereo information form the SMILES string?
    :type stereo: bool
    :rtype: str
    """
    if not stereo:
        smi = without_stereo(smi)

    smi = without_resonance_stereo(smi)
    rdm = rdkit_.from_smiles(smi)
    ich = rdkit_.to_inchi(rdm)
    return ich


def chi(smi):
    """Convert a SMILES string to an AMChI or InChI string.

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
    """Generate a molecular graph from a SMILES string.

    :param smi: SMILES string
    :type smi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param local_stereo: assign local stereo parities?
    :type local_stereo: bool
    :rtype: automol molecular graph
    """
    smis = split(smi)
    gras = [_connected_graph(s, stereo=stereo, local_stereo=local_stereo) for s in smis]
    gra = automol.graph.base.union_from_sequence(gras, shift_keys=True)
    return gra


def _connected_graph(smi, stereo=True, local_stereo=False):
    """Generate a connected molecular graph from a connected SMILES string.

    :param smi: SMILES string
    :type smi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param local_stereo: assign local stereo parities?
    :type local_stereo: bool
    :rtype: automol molecular graph
    """
    (
        symb_dct,
        bnd_ord_dct,
        atm_par_dct,
        bnd_par_dct,
    ) = parse_connected_molecule_properties(smi)
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
        sp2_bnd_keys = automol.graph.base.rigid_planar_bond_keys(gra)
        bnd_keys = ste_bnd_keys - sp2_bnd_keys
        gra = automol.graph.base.without_stereo(gra, bnd_keys=bnd_keys)

        if not local_stereo:
            # Convert from local to canonical stereo
            gra = automol.graph.base.from_local_stereo(gra)

    return gra


def geometry(smi, check=True):
    """Generate a molecular geometry from a SMILES string.

    :param smi: SMILES string
    :type smi: str
    :param check: check stereo and connectivity?
    :type check: bool
    :rtype: automol molecular geometry data structure
    """
    gra = graph(smi)
    geo = automol.graph.geometry(gra, check=check)
    return geo


def formula_string(smi):
    """Get the molecular formula string (Hill-sorted) from a SMILES string

    :param smi: SMILES string
    :type smi: str
    :rtype: str
    """
    ich = inchi(smi)
    return automol.inchi.base.formula_string(ich)


def recalculate_without_stereo(smi):
    """Recalculate a SMILES string, removing stereo if present

    :param smi: SMILES string
    :type smi: str
    :rtype: str
    """
    smi = without_stereo(smi)

    rdm = rdkit_.from_smiles(smi)
    smi = rdkit_.to_smiles(rdm)
    return smi


def rdkit_molecule(smi, stereo=True):
    """Convert a SMILES string to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(smi))

    :param smi: SMILES string
    :type smi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :returns: the RDKit molecule
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(smi, stereo=stereo)
    return automol.graph.rdkit_molecule(gra, stereo=stereo)


def svg_string(smi, stereo=True):
    """Convert a SMILES string into an SVG string for visualization

    :param smi: SMILES string
    :type smi: str
    :rtype: str
    """
    if not stereo:
        smi = without_stereo(smi)

    rdm = rdkit_.from_smiles(smi)
    svg_str = rdkit_.to_svg_string(rdm)
    return svg_str


def rdkit_reaction(rsmis, psmis, stereo=True, res_stereo=False):
    """Convert reactant and product graphs to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(pgras, rgras))

        :param rsmis: SMILES strings for the reactants
        :param psmis: SMILES strings for the products
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rgras = [graph(s, stereo=stereo) for s in rsmis]
    pgras = [graph(s, stereo=stereo) for s in psmis]
    return automol.graph.rdkit_reaction(
        rgras, pgras, stereo=stereo, res_stereo=res_stereo
    )


def display(smi, stereo=True):
    """Display graph to IPython using the RDKit visualizer

    :param smi: SMILES string
    :type smi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(smi, stereo=stereo)
    automol.graph.display(gra, stereo=stereo)


def display_reaction(rsmis, psmis, stereo=True):
    """Display reaction to IPython using the RDKit visualizer

    :param rsmis: SMILES strings for the reactants
    :param psmis: SMILES strings for the products
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    rgras = [graph(s, stereo=stereo) for s in rsmis]
    pgras = [graph(s, stereo=stereo) for s in psmis]
    automol.graph.display_reaction(rgras, pgras, stereo=stereo)


# helpers
def _compare(smi1, smi2):
    """Check if two SMILES strings are similar.

    :param smi1: SMILES string 1
    :type smi1: str
    :param smi2: SMILES string 2
    :type smi2: str
    :rtype: bool
    """
    return _canonicalize(smi1) == _canonicalize(smi2)


def _canonicalize(smi):
    """Convert a SMILES string into its canonical form.

    :param smi: SMILES string
    :type smi: str
    :rtype: str
    """
    return rdkit_.to_smiles(rdkit_.from_smiles(smi))
