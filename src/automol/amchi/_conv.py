""" Level 4 functions depending on other basic types (geom, graph)
"""

import itertools
import numbers

from .. import form, geom, inchi_key
from .. import graph as graph_
from ..extern import rdkit_
from ..util import dict_
from .base import (
    atom_stereo_parities,
    bond_stereo_parities,
    bonds,
    breaking_bond_keys,
    equivalent,
    forming_bond_keys,
    formula,
    formula_layer,
    has_stereo,
    hydrogen_valences,
    is_enantiomer,
    is_inchi,
    is_inverted_enantiomer,
    is_racemic,
    is_reversed_ts,
    isotope_layers,
    main_layers,
    racemic,
    split,
    standard_form,
    symbols,
    with_inchi_prefix,
)


# # conversions
def amchi_key(chi):
    """Generate a ChIKey from a ChI string.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    chi = with_inchi_prefix(chi)
    return rdkit_.inchi_to_inchi_key(chi)


def chemkin_name(chi: str, root_name: str | None = None) -> str:
    """Generate a CHEMKIN name from a ChI string.

    Todo: Simplified names for common species. Anything that can be described as a
    simple backbone with functional groups gets a simplified name. Examples:

        C3y1 = [CH2]CC = propan-1-yl
        C3p1 = [O]OCCC = propan-1-peroxy radical
        C3h1y2 = OOC[CH]C = propan-1-hydroperoxy-2-yl radical
        C3e1 = C=CC = prop-1-ene
        C3a1 = O=CCC = propan-1-al (always 1 for aldehydes)
        C3k2 = CC(=O)C = propan-2-one (never 1 for ketones)
        C4e12 = O1CC1CC = 1,2-epoxybutane
        C5e13 = O1CCC1CC = 1,3-epoxypentane

    :param chi: ChI string
    :param root_name: A root for the CHEMKIN name, to add stereo if needed
    :return: The CHEMKIN name
    """
    if root_name is not None:
        ste_str = stereo_digest(chi, max_len=16 - len(root_name))
        return root_name + ste_str

    return digest(chi, conn_max_len=10, ste_max_len=6)


def digest(chi: str, conn_max_len: int = 10, ste_max_len: int = 6) -> str:
    """Generate a short string summarizing connectivity information.

    :param chi: ChI string
    :param conn_length: Connectivity hash length, up to 14, defaults to 10
    :param stereo_length: Stereo hash length, up to 8, defaults to 5
    :return: The short hash
    """
    conn_str = connectivity_digest(chi, max_len=conn_max_len)
    ste_str = stereo_digest(chi, max_len=ste_max_len)
    return conn_str + ste_str


def connectivity_digest(chi: str, max_len: int = 10) -> str:
    """Generate a short string summarizing connectivity information.

    :param chi: ChI string
    :param max_len: The maximum allowed length of the digest
    :return: The short hash
    """
    # hash length by heavy atom count
    hash_len_dct = {0: 0, 1: 0, 2: 4, 3: 5, 4: 5, 5: 6, 6: 6, 7: 6, 8: 6}

    # Determine appropriate hash length
    fml = formula(chi)
    nheavy = form.heavy_atom_count(fml)
    hash_len = hash_len_dct[nheavy]

    # Form the hash string
    chk = amchi_key(chi)
    hash_str = inchi_key.first_hash(chk)[:hash_len].lower()

    # Determine approprriate formula length
    form_len = max_len - hash_len if max_len > hash_len else 0

    # Form the formula string
    form_strs = [form.string(fml), form.string(fml, hyd=False), ""]
    form_str = next(s for s in form_strs if len(s) <= form_len)

    return form_str + hash_str


def stereo_digest(chi: str, racem: bool = False, max_len: int | None = None) -> str:
    """Generate a short string summarizing stereo information.

    :param name: The CHEMKIN name
    :param chi: The ChI string
    :param racem: Use a racemic suffix?, defaults to False
    :param max_len: The maximum allowed length of the digest
    :return: The CHEMKIN name, with stereo suffix
    """
    if not has_stereo(chi):
        return ""

    # Build the stereo digest
    bpar_dct = bond_stereo_parities(chi, ordered_key=True)
    apar_dct = atom_stereo_parities(chi)
    is_inv = is_inverted_enantiomer(chi)
    is_enant = is_enantiomer(chi)
    bpars = list(map(bpar_dct.get, sorted(bpar_dct)))
    apars = list(map(apar_dct.get, sorted(apar_dct)))
    bnd_str = "".join("e" if p else "z" for p in bpars)
    atm_str = "".join("s" if p else "r" for p in apars)
    inv_str = ""
    if is_enant:
        racem = True if is_racemic(chi) else racem
        inv_str = "R" if racem else ("1" if is_inv else "0")
    ste_str = "".join((bnd_str, atm_str, inv_str))

    # If longer than the max length, use a hash instead
    if max_len is not None and len(bnd_str + atm_str) > max_len - 1:
        hash2 = inchi_key.second_hash(amchi_key(racemic(chi)))
        ste_str = hash2[: max_len - 1] + inv_str

    return ste_str


def chi_(chi: str) -> str:
    """Convert AMChI string to InChI, if appropriate.

    :param chi: AMChI string
    :return: An InChI or AMChI string encoding the same species
    """
    gra = graph(chi, local_stereo=True)
    ich = graph_.inchi(gra, stereo=True, local_stereo=True)
    if not graph_.inchi_is_bad(gra, ich):
        return ich
    return chi


def smiles(chi, res_stereo=True):
    """Convert a ChI string into a SMILES string.

    :param chi: ChI string
    :type chi: str
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: the SMILES string
    :rtype: str
    """
    chis = split(chi)
    smis = [_connected_smiles(c, res_stereo=res_stereo) for c in chis]
    smi = ".".join(smis)
    return smi


def _connected_smiles(chi, res_stereo=True):
    """Convert a single-component ChI string into a SMILES string.

    :param chi: ChI string
    :type chi: str
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: the SMILES string
    :rtype: str
    """
    gra = _connected_graph(chi, stereo=True, local_stereo=True)
    smi = graph_.smiles(gra, stereo=True, local_stereo=True, res_stereo=res_stereo)
    return smi


def graph(chi, stereo=True, local_stereo=False):
    """Generate a molecular graph from a ChI string.

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param local_stereo: assign local stereo parities?
    :type local_stereo: bool
    :rtype: automol molecular graph
    """
    chis = split(chi)
    gras = [_connected_graph(c, stereo=stereo, local_stereo=local_stereo) for c in chis]
    gra = graph_.union_from_sequence(gras, shift_keys=True)
    if stereo:
        gra = graph_.with_explicit_stereo_hydrogens(gra)
    return gra


def _connected_graph(chi, stereo=True, local_stereo=False):
    """Generate a connected molecular graph from a single-component ChI string

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param local_stereo: assign local stereo parities?
    :type local_stereo: bool
    :rtype: automol molecular graph
    """
    symb_dct = symbols(chi)
    bnd_keys = bonds(chi)
    atm_imp_hyd_dct = hydrogen_valences(chi)

    if isotope_layers(chi):
        raise NotImplementedError("Isotopic graph conversion not implemented")

    if stereo:
        bnd_ste_par_dct = bond_stereo_parities(chi)
        atm_ste_par_dct = atom_stereo_parities(chi)
    else:
        atm_ste_par_dct = None
        bnd_ste_par_dct = None

    is_inv = is_inverted_enantiomer(chi)

    # TS graph
    brk_bkeys = breaking_bond_keys(chi)
    frm_bkeys = forming_bond_keys(chi)
    is_rev = is_reversed_ts(chi)

    all_keys = set(symb_dct.keys()) | set(itertools.chain(*bnd_keys))
    symb_dct = dict_.by_key(symb_dct, all_keys, fill_val="H")
    gra = graph_.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
    )

    if brk_bkeys or frm_bkeys:
        gra = graph_.ts_graph(gra, frm_bnd_keys=frm_bkeys, brk_bnd_keys=brk_bkeys)

    if is_inv is True:
        gra = graph_.invert_atom_stereo_parities(gra)

    # Note: If this is an InChI string, local != canonical stereo!
    if has_stereo(chi) and not local_stereo and (is_inv or is_inchi(chi)):
        gra = graph_.from_local_stereo(gra)

    if is_rev is True:
        gra = graph_.ts_reverse(gra)

    return gra


def geometry(chi, check=True, log=False):
    """Generate a molecular geometry from a ChI string.

    :param chi: ChI string
    :type chi: str
    :param check: check stereo and connectivity?
    :type check: bool
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :rtype: automol molecular geometry data structure
    """
    gra = graph(chi)

    if log:
        print("Graph generated from {chi}:")
        print(gra)

    geo = graph_.geometry(gra, check=check, log=log)
    return geo


def zmatrix(chi, check=True):
    """Generate a z-matrix from an ChI string.

    :param chi: ChI string
    :type chi: str
    :param check: check stereo and connectivity?
    :type check: bool
    :rtype: automol z-matrix data structure
    """

    geo = geometry(chi, check=check)
    zma = geom.zmatrix(geo)
    return zma


def rdkit_molecule(chi, stereo=True):
    """Convert a ChI string to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(chi))

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :returns: the RDKit molecule
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(chi, stereo=stereo)
    return graph_.rdkit_molecule(gra, stereo=stereo)


def rdkit_reaction(rchis, pchis, stereo=True, res_stereo=False):
    """Convert reactant and product graphs to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(pgras, rgras))

        :param rchis: ChI strings for the reactants
        :param pchis: ChI strings for the products
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rgras = [graph(s, stereo=stereo) for s in rchis]
    pgras = [graph(s, stereo=stereo) for s in pchis]
    return graph_.rdkit_reaction(rgras, pgras, stereo=stereo, res_stereo=res_stereo)


def display(chi, stereo=True, exp=False, label=False, label_dct=None):
    """Display molecule to IPython using the RDKit visualizer

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(chi, stereo=stereo)
    graph_.display(gra, stereo=stereo, exp=exp, label=label, label_dct=label_dct)


def display_reaction(rchis, pchis, stereo=True):
    """Display reaction to IPython using the RDKit visualizer

    :param rchis: ChI strings for the reactants
    :param pchis: ChI strings for the products
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    rgras = [graph(s, stereo=stereo) for s in rchis]
    pgras = [graph(s, stereo=stereo) for s in pchis]
    graph_.display_reaction(rgras, pgras, stereo=stereo)


# # derived properties
def is_complete(chi):
    """Determine if the ChI string is complete
    (has all stereo-centers assigned).

    Currently only checks species that does not have any
    resonance structures.

    :param chi: ChI string
    :type chi: str
    :rtype: bool
    """
    gra = graph(chi, stereo=False)
    is_missing_stereo = bool(graph_.unassigned_stereocenter_keys(gra))

    return equivalent(chi, standard_form(chi)) and not (
        has_stereo(chi) ^ is_missing_stereo
    )


def is_valid_multiplicity(chi, mul):
    """is this multiplicity compatible with this amchi string?

    :param chi: ChI string
    :type chi: str
    :param mul: multiplicity
    :type mul: int
    :returns: validity of amchi multiplicity
    :rtype: bool
    """
    assert isinstance(mul, numbers.Integral)
    return mul in graph_.possible_spin_multiplicities(graph(chi, stereo=False))


def guess_spin(chi: str) -> int:
    """Guess the spin of a ChI string.

    Apart from a special list of common molecules, simply guesses 0 or 1 depending on
    whether the eletron count is even or odd, respectively.

    Could be made smarter by converting to a graph, but that would be more expensive.

    :param chi: A ChI string
    :return: The ground-state spin guess
    """
    # First, check if this is a case with a hardcoded value
    lookup_dct = {
        ("O",): 2,
        ("O2", ("c", "1-2")): 2,
    }
    key = (formula_layer(chi), *sorted(main_layers(chi).items()))
    val = lookup_dct.get(key, None)
    if val is not None:
        return val

    # Otherwise, guess based on the formula
    fml = formula(chi)
    return form.electron_count(fml) % 2


# # derived transformations
def add_stereo(chi):
    """Add stereochemistry to a ChI string converting to/from geometry.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    geo = geometry(chi)
    chi = geom.amchi(geo, stereo=True)
    return chi


def expand_stereo(chi: str, enant: bool = True, strained: bool = False) -> list[str]:
    """Obtain all possible stereoisomers of a ChI string.

    :param chi: ChI string
    :param enant: Include all enantiomers, or only canonical ones?
    :param strained: Include strained stereoisomers?
    """
    gra = graph(chi, stereo=False)
    sgrs = graph_.expand_stereo(gra, enant=enant, strained=strained, symeq=False)
    ste_chis = [graph_.amchi(sgr, stereo=True) for sgr in sgrs]
    return ste_chis
