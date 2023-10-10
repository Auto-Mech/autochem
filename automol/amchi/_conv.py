""" Level 4 functions depending on other basic types (geom, graph)
"""
import automol.graph
import automol.geom
from automol import error
from automol.extern import rdkit_
from automol.amchi.base import isotope_layers
from automol.amchi.base import symbols
from automol.amchi.base import bonds
from automol.amchi.base import hydrogen_valences
from automol.amchi.base import atom_stereo_parities
from automol.amchi.base import bond_stereo_parities
from automol.amchi.base import is_inverted_enantiomer
from automol.amchi.base import has_stereo
from automol.amchi.base import split
from automol.amchi.base import with_inchi_prefix
from automol.amchi.base import equivalent
from automol.amchi.base import standard_form


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
    smi = automol.graph.smiles(
        gra, stereo=True, local_stereo=True, res_stereo=res_stereo
    )
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
    gra = automol.graph.union_from_sequence(gras, shift_keys=True)
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

    gra = automol.graph.from_data(
        atm_symb_dct=symb_dct,
        bnd_keys=bnd_keys,
        atm_imp_hyd_dct=atm_imp_hyd_dct,
        atm_ste_par_dct=atm_ste_par_dct,
        bnd_ste_par_dct=bnd_ste_par_dct,
    )

    if is_inv is True:
        gra = automol.graph.reflect_local_stereo(gra)

    if has_stereo(chi) and not local_stereo:
        gra = automol.graph.from_local_stereo(gra)

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

    geo = automol.graph.geometry(gra, check=check, log=log)
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
    gra = automol.graph.explicit(gra)

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
                gra_ = automol.geom.graph(geo, stereo=False)
                if automol.graph.isomorphism(gra, gra_, stereo=False):
                    ret_geos.append(geo)
            else:
                # There is stereo.
                # First, check connectivity.
                gra_ = automol.geom.graph(geo)
                geo_idx_dct = automol.graph.isomorphism(gra, gra_, stereo=False)

                if geo_idx_dct is not None:
                    # Enforce correct stereo parities. This is necessary for
                    # resonance bond stereo.
                    geo = automol.graph.stereo_corrected_geometry(
                        gra, geo, geo_idx_dct=geo_idx_dct, local_stereo=True
                    )

                    # Check if the assignment worked.
                    gra_ = automol.graph.set_stereo_from_geometry(gra_, geo)
                    if automol.graph.isomorphism(gra, gra_):
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
    zma = automol.geom.zmatrix(geo)
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
    return automol.graph.rdkit_molecule(gra, stereo=stereo)


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
    return automol.graph.rdkit_reaction(
        rgras, pgras, stereo=stereo, res_stereo=res_stereo
    )


def display(chi, stereo=True):
    """Display molecule to IPython using the RDKit visualizer

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    gra = graph(chi, stereo=stereo)
    automol.graph.display(gra, stereo=stereo)


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
    automol.graph.display_reaction(rgras, pgras, stereo=stereo)


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
    ste_atm_keys = automol.graph.stereogenic_atom_keys(gra)
    ste_bnd_keys = automol.graph.stereogenic_bond_keys(gra)
    graph_has_stereo = bool(ste_atm_keys or ste_bnd_keys)

    _complete = equivalent(chi, standard_form(chi)) and not (
        has_stereo(chi) ^ graph_has_stereo
    )

    return _complete


# # derived transformations
def add_stereo(chi):
    """Add stereochemistry to a ChI string converting to/from geometry.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    geo = geometry(chi)
    chi = automol.geom.amchi(geo, stereo=True)
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
    sgrs = automol.graph.expand_stereo(gra, enant=enant, symeq=False)
    ste_chis = [automol.graph.amchi(sgr, stereo=True) for sgr in sgrs]
    return ste_chis
