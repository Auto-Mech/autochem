""" Level 4 functions depending on other basic types (geom, graph)
"""
import itertools
import numbers

from automol import error, geom, graph as graph_
from automol.amchi.base import (
    atom_stereo_parities,
    bond_stereo_parities,
    bonds,
    breaking_bond_keys,
    equivalent,
    forming_bond_keys,
    has_stereo,
    hydrogen_valences,
    is_inchi,
    is_inverted_enantiomer,
    is_reversed_ts,
    isotope_layers,
    split,
    standard_form,
    symbols,
    with_inchi_prefix,
)
from automol.extern import rdkit_
from automol.util import dict_


# # conversions
def amchi_key(chi):
    """Generate a ChIKey from a ChI string.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    chi = with_inchi_prefix(chi)
    return rdkit_.inchi_to_inchi_key(chi)


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


def conformers(chi, nconfs=1, check=True, accept_fewer=False):
    """Generate a connected molecular geometry from a single-component ChI
    string.

    :param chi: ChI string
    :type chi: str
    :param nconfs: number of conformer structures to generate
    :type nconfs: int
    :param check: check stereo and connectivity?
    :type check: bool
    :param accept_fewer: accept fewer than nconfs conformers?
    :type accept_fewer: bool
    :rtype: automol molecular geometry
    """
    # Convert graph to local stereo to avoid multiple recanonicalizations
    gra = _connected_graph(chi, stereo=True, local_stereo=True)
    gra = graph_.explicit(gra)

    smi = _connected_smiles(chi, res_stereo=False)
    has_ste = has_stereo(chi)

    try:
        rdm = rdkit_.from_smiles(smi)
        geos = rdkit_.to_conformers(rdm, nconfs=nconfs)
    except (RuntimeError, TypeError, ValueError) as err:
        raise error.FailedGeometryGenerationError("Failed AMChI:", chi) from err

    ret_geos = []

    # If the ChI has stereo, enforce correct stereo on the geometry.
    if check:
        for geo in geos:
            if not has_ste:
                # If there wasn't stereo, only check connectivity
                gra_ = geom.graph(geo, stereo=False)
                if graph_.isomorphism(gra, gra_, stereo=False):
                    ret_geos.append(geo)
            else:
                # There is stereo.
                # First, check connectivity.
                gra_ = geom.graph(geo)
                geo_idx_dct = graph_.isomorphism(gra, gra_, stereo=False)

                if geo_idx_dct is not None:
                    # Enforce correct stereo parities. This is necessary for
                    # resonance bond stereo.
                    geo = graph_.stereo_corrected_geometry(
                        gra, geo, geo_idx_dct=geo_idx_dct, local_stereo=True
                    )

                    # Check if the assignment worked.
                    gra_ = graph_.set_stereo_from_geometry(gra_, geo)
                    if graph_.isomorphism(gra, gra_):
                        ret_geos.append(geo)

    if len(ret_geos) < nconfs and not accept_fewer:
        raise error.FailedGeometryGenerationError("Failed AMChI:", chi)

    return ret_geos


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


def display(chi, stereo=True):
    """Display molecule to IPython using the RDKit visualizer

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(chi, stereo=stereo)
    graph_.display(gra, stereo=stereo)


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


def expand_stereo(chi, enant=True):
    """Obtain all possible stereoisomers of a ChI string.

    :param chi: ChI string
    :type chi: str
    :param enant: Include all enantiomers, or only canonical ones?
    :type enant: bool
    :rtype: list[str]
    """
    gra = graph(chi, stereo=False)
    sgrs = graph_.expand_stereo(gra, enant=enant, symeq=False)
    ste_chis = [graph_.amchi(sgr, stereo=True) for sgr in sgrs]
    return ste_chis
