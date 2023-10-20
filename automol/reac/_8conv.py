""" Conversion functions
"""
from typing import List

import IPython

from automol import amchi, geom, graph
from automol import chi as chi_
from automol import smiles as smiles_
from automol import zmat as zmat_
from automol.extern import rdkit_
from automol.reac._0core import (
    Reaction,
    apply_zmatrix_conversion,
    product_graphs,
    reactant_graphs,
    set_structures,
    ts_graph,
    without_stereo,
)
from automol.reac._2stereo import reflect
from automol.reac._3find import find
from automol.reac._4struc import ts_geometry_from_reactants, with_structures
from automol.reac._5zmat import ts_zmatrix


# Conversion stuff
def amchis(rxn: Reaction, stereo=True):
    """Convert the reaction object to AMChIs

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: AMChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_chis = tuple(graph.amchi(gra, stereo=stereo) for gra in reactant_graphs(rxn))
    prd_chis = tuple(graph.amchi(gra, stereo=stereo) for gra in product_graphs(rxn))
    return (rct_chis, prd_chis)


def ts_amchi(rxn: Reaction, stereo=True) -> str:
    """Get the AMChI for the reaction TS

    :param rxn: The reaction object
    :type rxn: Reaction
    :param stereo: Include stereo?
    :type stereo: bool
    :returns: AMChI string for the TS
    :rtype: str
    """
    if not stereo:
        rxn = without_stereo(rxn)

    tsg = ts_graph(rxn)
    ts_chi = graph.amchi(tsg)
    return ts_chi


def inchis(rxn: Reaction, stereo=True):
    """Convert the reaction object to InChIs

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: InChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_ichs = tuple(graph.inchi(gra, stereo=stereo) for gra in reactant_graphs(rxn))
    prd_ichs = tuple(graph.inchi(gra, stereo=stereo) for gra in product_graphs(rxn))
    return (rct_ichs, prd_ichs)


def chis(rxn: Reaction, stereo=True):
    """Convert the reaction object to ChIs

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: ChI strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_chis = tuple(graph.chi(gra, stereo=stereo) for gra in reactant_graphs(rxn))
    prd_chis = tuple(graph.chi(gra, stereo=stereo) for gra in product_graphs(rxn))
    return (rct_chis, prd_chis)


def smiles(rxn: Reaction, stereo=True, res_stereo=True, exp_singles=False):
    """Convert the reaction object to SMILESs

    :param rxn: the reaction object
    :param stereo: Include stereo?
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :param exp_singles: Use explicit '-' for single bonds?
    :type exp_singles: bool
    :returns: SMILES strings for the reactants and products
    :rtype: (tuple[str], tuple[str])
    """
    rct_smis = tuple(
        graph.smiles(gra, stereo=stereo, res_stereo=res_stereo, exp_singles=exp_singles)
        for gra in reactant_graphs(rxn)
    )
    prd_smis = tuple(
        graph.smiles(gra, stereo=stereo, res_stereo=res_stereo, exp_singles=exp_singles)
        for gra in product_graphs(rxn)
    )
    return (rct_smis, prd_smis)


def reaction_smiles(rxn) -> str:
    """Convert the Reaction object to a reaction SMILES string

    :param rxn: The reaction object
    :type rxn: Reaction
    :returns: The reaction SMILES
    :rtype: str
    """
    rct_smis, prd_smis = smiles(rxn)
    rxn_smi = smiles_.reaction(rct_smis, prd_smis)
    return rxn_smi


def rdkit_reaction(rxn: Reaction, stereo=True, res_stereo=False):
    """Convert a reaction object to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(gra))

        :param rxn: the reaction object
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :param rxn: the reaction object
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rct_smis, prd_smis = smiles(
        rxn, stereo=stereo, res_stereo=res_stereo, exp_singles=True
    )
    rcts_smi = ".".join(rct_smis)
    prds_smi = ".".join(prd_smis)
    rxn_smi = f"{rcts_smi}>>{prds_smi}"
    return rdkit_.from_smarts(rxn_smi)


def display(rxn: Reaction):
    """Display reaction object to IPython using the RDKit visualizer

    :param rxn: the reaction object
    :returns: None
    """
    rdkit_.turn_3d_visualization_off()
    return IPython.display.display(rdkit_reaction(rxn))


# Check if the reaction is canonical
def is_canonical_enantiomer(srxn: Reaction):
    """Does this reaction have a canonical combination of enantiomers?

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    rct_chis, prd_chis = amchis(srxn)
    return amchi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)


def canonical_enantiomer(srxn: Reaction):
    """Convert this reaction into a canonical combination of enantiomers

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    if not is_canonical_enantiomer(srxn):
        srxn = reflect(srxn)
    return srxn


# Get a reaction object with structures from various identifiers
def from_chis(
    rct_chis, prd_chis, stereo=False, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from ChIs

    :param rct_chis: The reactant ChI (InChI or AMChI) strings
    :type rct_chis: list[str]
    :param prd_chis: The product ChI (InChI or AMChI) strings
    :type prd_chis: list[str]
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to None
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(graph.explicit, map(chi_.graph, rct_chis)))
    prd_gras = tuple(map(graph.explicit, map(chi_.graph, prd_chis)))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_smiles(
    rct_smis, prd_smis, stereo=False, struc_typ: str = None
) -> List[Reaction]:
    """Get reaction objects from SMILES

    :param rct_smis: The reactant SMILES strings
    :type rct_smis: list[str]
    :param prd_smis: The product SMILES strings
    :type prd_smis: list[str]
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to None
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(graph.explicit, map(smiles_.graph, rct_smis)))
    prd_gras = tuple(map(graph.explicit, map(smiles_.graph, prd_smis)))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_zmatrices(
    rct_zmas, prd_zmas, stereo=False, struc_typ: str = "zmat"
) -> List[Reaction]:
    """Get reaction objects from z-matrices

    :param rct_zmas: The reactant z-matrices
    :type rct_zmas: list of automol z-matrices
    :param prd_zmas: The product z-matrices
    :type prd_zmas: list of automol z-matrices
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to "zmat"
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(zmat_.graph, rct_zmas))
    prd_gras = tuple(map(zmat_.graph, prd_zmas))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def from_geometries(
    rct_geos, prd_geos, stereo=False, struc_typ: str = "geom"
) -> List[Reaction]:
    """Get reaction objects from geometries

    :param rct_geos: The reactant geometries
    :type rct_geos: list of automol geometries
    :param prd_geos: The product geometries
    :type prd_geos: list of automol geometries
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :param struc_typ: Add structures of this type; defaults to "geom"
    :type struc_typ: str, optional
    :returns: A series of reaction objects
    :rtype: List[Reaction]
    """
    rct_gras = tuple(map(geom.graph, rct_geos))
    prd_gras = tuple(map(geom.graph, prd_geos))
    rxns = find(rct_gras, prd_gras, stereo=stereo)

    if struc_typ is not None:
        rxns = tuple(with_structures(r, struc_typ) for r in rxns)

    return rxns


def with_structures_from_chi(rct_chis, prd_chis, zmat=False, stereo=False):
    """Get reaction objects with geometry/z-matrix structures from ChIs

    :param rct_chis: The reactant ChI (InChI or AMChI) strings
    :type rct_chis: list[str]
    :param prd_chis: The product ChI (InChI or AMChI) strings
    :type prd_chis: list[str]
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    rct_geos = list(map(chi_.geometry, rct_chis))
    prd_geos = list(map(chi_.geometry, prd_chis))

    return with_structures_from_geometry(rct_geos, prd_geos, zmat=zmat, stereo=stereo)


def with_structures_from_smiles(rct_smis, prd_smis, zmat=False, stereo=False):
    """Get reaction objects with geometry/z-matrix structures from SMILES

    :param rct_smis: The reactant SMILES strings
    :type rct_smis: list[str]
    :param prd_smis: The product SMILES strings
    :type prd_smis: list[str]
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    # rct_geos = list(map(smiles_.geometry, rct_smis))
    # prd_geos = list(map(smiles_.geometry, prd_smis))

    rct_chis = list(map(smiles_.chi, rct_smis))
    prd_chis = list(map(smiles_.chi, prd_smis))
    rct_geos = list(map(chi_.geometry, rct_chis))
    prd_geos = list(map(chi_.geometry, prd_chis))

    return with_structures_from_geometry(rct_geos, prd_geos, zmat=zmat, stereo=stereo)


def with_structures_from_zmatrix(rct_zmas, prd_zmas, zmat=False, stereo=False):
    """Get reaction objects with geometry/z-matrix structures from SMILES

    :param rct_zmas: The reactant z-matrices
    :type rct_zmas: list of automol z-matrices
    :param prd_zmas: The product z-matrices
    :type prd_zmas: list of automol z-matrices
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :returns: A series of tuples containing, in order, the reaction object, the
        TS structure, the (sorted) reactant structures, and the (sorted)
        product structures. The latter three are returned as either automol
        geometry objects or, if `zmatrix` was set to `True, automol z-matrix
        objects.
    """

    rct_geos = list(map(zmat_.geometry, rct_zmas))
    prd_geos = list(map(zmat_.geometry, prd_zmas))

    return with_structures_from_geometry(rct_geos, prd_geos, zmat=zmat, stereo=stereo)


def with_structures_from_geometry(
    rct_geos, prd_geos, zmat=False, stereo=False
) -> List[Reaction]:
    """Get reaction objects with geometry/z-matrix structures from geometries

    :param rct_geos: The reactant geometries
    :type rct_geos: list of automol geometries
    :param prd_geos: The product geometries
    :type prd_geos: list of automol geometries
    :param zmat: Return z-matrix structures instead of cartesian geometries?
    :type zmat: bool
    :param stereo: Include stereoassignments?
    :type stereo: bool
    :returns: A series of Reaction objects containing structures
    :rtype: List[Reaction]
    """

    # Identify the reaction based on the reactants and products
    rct_gras = list(map(geom.graph, rct_geos))
    prd_gras = list(map(geom.graph, prd_geos))

    if not stereo:
        rct_gras = list(map(graph.without_stereo, rct_gras))
        prd_gras = list(map(graph.without_stereo, prd_gras))

    rxns = find(rct_gras, prd_gras)

    # Obtain the reaction objects and structures to return
    ret = []
    for rxn in rxns:
        ts_geo = ts_geometry_from_reactants(rxn, rct_geos, log=False)
        rxn = set_structures(rxn, ts_geo, rct_geos, prd_geos)
        # Determine which geometries to store
        if not zmat:
            ret.append(rxn)
        else:
            ts_zma, dc_ = ts_zmatrix(rxn, ts_geo)
            zrxn = apply_zmatrix_conversion(rxn, dc_)
            rct_zmas = tuple(map(geom.zmatrix, rct_geos))
            prd_zmas = tuple(map(geom.zmatrix, prd_geos))

            ret += ((zrxn, ts_zma, rct_zmas, prd_zmas),)

    # Set to None if no objects found
    if not ret:
        ret = None

    return ret
